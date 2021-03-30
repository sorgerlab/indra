__all__ = ['Agent', 'BoundCondition', 'MutCondition', 'ModCondition',
           'ActivityCondition', 'default_ns_order']


import logging
from collections import OrderedDict as _o
from indra.statements.statements import modtype_conditions, modtype_to_modclass
from indra.statements.validate import assert_valid_db_refs
from .concept import Concept
from .resources import get_valid_residue, activity_types, amino_acids


logger = logging.getLogger(__name__)


default_ns_order = ['FPLX', 'UPPRO', 'HGNC', 'UP', 'CHEBI', 'GO', 'MESH',
                    'MIRBASE', 'DOID', 'HP', 'EFO']


class Agent(Concept):
    """A molecular entity, e.g., a protein.

    Parameters
    ----------
    name : str
        The name of the agent, preferably a canonicalized name such as an
        HGNC gene name.
    mods : list of :py:class:`ModCondition`
        Modification state of the agent.
    bound_conditions : list of :py:class:`BoundCondition`
        Other agents bound to the agent in this context.
    mutations : list of :py:class:`MutCondition`
        Amino acid mutations of the agent.
    activity : :py:class:`ActivityCondition`
        Activity of the agent.
    location : str
        Cellular location of the agent. Must be a valid name (e.g. "nucleus")
        or identifier (e.g. "GO:0005634")for a GO cellular compartment.
    db_refs : dict
        Dictionary of database identifiers associated with this agent.
    """
    def __init__(self, name, mods=None, activity=None,
                 bound_conditions=None, mutations=None,
                 location=None, db_refs=None):
        super(Agent, self).__init__(name, db_refs=db_refs)

        if mods is None:
            self.mods = []
        # Promote to list
        elif isinstance(mods, ModCondition):
            self.mods = [mods]
        else:
            self.mods = mods

        if bound_conditions is None:
            self.bound_conditions = []
        # Promote to list
        elif isinstance(bound_conditions, BoundCondition):
            self.bound_conditions = [bound_conditions]
        else:
            self.bound_conditions = bound_conditions

        if mutations is None:
            self.mutations = []
        elif isinstance(mutations, MutCondition):
            self.mutations = [mutations]
        else:
            self.mutations = mutations

        self.activity = activity
        self.location = location

    @classmethod
    def from_refs(cls, name, db_refs, **kwargs) -> 'Agent':
        """Create an agent from db_refs."""
        from ..ontology.standardize import standardize_name_db_refs
        standard_name, db_refs = standardize_name_db_refs(db_refs)
        if standard_name:
            name = standard_name
        assert_valid_db_refs(db_refs)
        return cls(name, db_refs=db_refs, **kwargs)

    def matches_key(self):
        """Return a key to identify the identity and state of the Agent."""
        key = (self.entity_matches_key(),
               self.state_matches_key())
        return str(key)

    def entity_matches_key(self):
        """Return a key to identify the identity of the Agent not its state.

        The key is based on the preferred grounding for the Agent, or if not
        available, the name of the Agent is used.

        Returns
        -------
        str
            The key used to identify the Agent.
        """
        db_ns, db_id = self.get_grounding()
        if db_ns and db_id:
            return str((db_ns, db_id))
        return self.name

    def state_matches_key(self):
        """Return a key to identify the state of the Agent."""
        # NOTE: Making a set of the mod matches_keys might break if
        # you have an agent with two phosphorylations at serine
        # with unknown sites.
        act_key = (self.activity.matches_key() if self.activity else None)
        key = (sorted([m.matches_key() for m in self.mods]),
               sorted([m.matches_key() for m in self.mutations]),
               act_key, self.location,
               len(self.bound_conditions),
               tuple((bc.agent.matches_key(), bc.is_bound)
                     for bc in sorted(self.bound_conditions,
                                      key=lambda x: x.agent.name)))
        return str(key)

    # Function to get the namespace to look in
    def get_grounding(self, ns_order=None):
        """Return a tuple of a preferred grounding namespace and ID.

        Returns
        -------
        tuple
            A tuple whose first element is a grounding namespace (HGNC,
            CHEBI, etc.) and the second element is an identifier in the
            namespace. If no preferred grounding is available, a tuple of
            Nones is returned.
        """
        return get_grounding(self.db_refs, ns_order=ns_order)

    def isa(self, other, ontology):
        # Get the namespaces for the comparison
        (self_ns, self_id) = self.get_grounding()
        (other_ns, other_id) = other.get_grounding()
        # If one of the agents isn't grounded to a relevant namespace,
        # there can't be an isa relationship
        if not all((self_ns, self_id, other_ns, other_id)):
            return False
        # Check for isa relationship
        return ontology.isa_or_partof(self_ns, self_id, other_ns,
                                      other_id)

    def refinement_of(self, other, ontology, entities_refined=False):
        from indra.databases import go_client
        # Make sure the Agent types match
        if type(self) != type(other):
            return False

        # ENTITIES
        # Check that the basic entity of the agent either matches or is related
        # to the entity of the other agent. If not, no match.

        # If the entities, match, then we can continue
        if not (entities_refined or
                (self.entity_matches(other) or self.isa(other, ontology))):
            return False

        # BOUND CONDITIONS
        # Now check the bound conditions. For self to be a refinement of
        # other in terms of the bound conditions, it has to include all of the
        # bound conditions in the other agent, and add additional context.
        # TODO: For now, we do not check the bound conditions of the bound
        # conditions.
        # Iterate over the bound conditions in the other agent, and make sure
        # they are all matched in self.
        used_idx = set()
        for bc_other in other.bound_conditions:
            # Iterate over the bound conditions in self to find a match
            bc_found = False
            for idx, bc_self in enumerate(self.bound_conditions):
                if (idx not in used_idx) and \
                        (bc_self.is_bound == bc_other.is_bound) and \
                        bc_self.agent.refinement_of(bc_other.agent, ontology):
                    bc_found = True
                    used_idx.add(idx)
                    break
            # If we didn't find a match for this bound condition in self, then
            # no refinement
            if not bc_found:
                return False

        # MODIFICATIONS
        # Similar to the above, we check that self has all of the modifications
        # of other.
        # Here we need to make sure that a mod in self.mods is only matched
        # once to a mod in other.mods. Otherwise ('phoshporylation') would be
        # considered a refinement of ('phosphorylation', 'phosphorylation')
        matched_indices = []
        # This outer loop checks that each modification in the other Agent
        # is matched.
        for other_mod in other.mods:
            mod_found = False
            # We need to keep track of indices for this Agent's modifications
            # to make sure that each one is used at most once to match
            # the modification of one of the other Agent's modifications.
            for ix, self_mod in enumerate(self.mods):
                if self_mod.refinement_of(other_mod, ontology):
                    # If this modification hasn't been used for matching yet
                    if ix not in matched_indices:
                        # Set the index as used
                        matched_indices.append(ix)
                        mod_found = True
                        break
            # If we didn't find an exact match for this mod in other, then
            # no refinement
            if not mod_found:
                return False

        # MUTATIONS
        # Similar to the above, we check that self has all of the mutations
        # of other.
        matched_indices = []
        # This outer loop checks that each mutation in the other Agent
        # is matched.
        for other_mut in other.mutations:
            mut_found = False
            # We need to keep track of indices for this Agent's mutations
            # to make sure that each one is used at most once to match
            # the mutation of one of the other Agent's mutations.
            for ix, self_mut in enumerate(self.mutations):
                if self_mut.refinement_of(other_mut):
                    # If this mutation hasn't been used for matching yet
                    if ix not in matched_indices:
                        # Set the index as used
                        matched_indices.append(ix)
                        mut_found = True
                        break
            # If we didn't find an exact match for this mut in other, then
            # no refinement
            if not mut_found:
                return False

        # LOCATION
        # If the other location is specified and this one is not then self
        # cannot be a refinement
        if self.location is None:
            if other.location is not None:
                return False
        # If both this location and the other one is specified, we check the
        # hierarchy.
        elif other.location is not None:
            # If the other location is part of this location then
            # self.location is not a refinement
            sl = go_client.get_go_id_from_label(self.location)
            ol = go_client.get_go_id_from_label(other.location)
            if not ontology.isa_or_partof('GO', sl, 'GO', ol):
                return False

        # ACTIVITY
        if self.activity is None:
            if other.activity is not None:
                return False
        elif other.activity is not None:
            if not self.activity.refinement_of(other.activity, ontology):
                return False

        # Everything checks out
        return True

    def equals(self, other):
        matches = (self.name == other.name) and \
                  (self.activity == other.activity) and \
                  (self.location == other.location) and \
                  (self.db_refs == other.db_refs)
        if len(self.mods) == len(other.mods):
            for s, o in zip(self.mods, other.mods):
                matches = matches and s.equals(o)
        else:
            return False
        if len(self.mutations) == len(other.mutations):
            for s, o in zip(self.mutations, other.mutations):
                matches = matches and s.equals(o)
        else:
            return False
        if len(self.bound_conditions) == len(other.bound_conditions):
            for s, o in zip(self.bound_conditions, other.bound_conditions):
                matches = matches and s.agent.equals(o.agent) and \
                          s.is_bound == o.is_bound
        else:
            return False

        return matches

    def to_json(self):
        json_dict = _o({'name': self.name})
        if self.mods:
            json_dict['mods'] = [mc.to_json() for mc in self.mods]
        if self.mutations:
            json_dict['mutations'] = [mc.to_json() for mc in self.mutations]
        if self.bound_conditions:
            json_dict['bound_conditions'] = [bc.to_json() for bc in
                                             self.bound_conditions]
        if self.activity is not None:
            json_dict['activity'] = self.activity.to_json()
        if self.location is not None:
            json_dict['location'] = self.location
        json_dict['db_refs'] = self.db_refs
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        name = json_dict.get('name')
        db_refs = json_dict.get('db_refs', {})
        mods = json_dict.get('mods', [])
        mutations = json_dict.get('mutations', [])
        activity = json_dict.get('activity')
        bound_conditions = json_dict.get('bound_conditions', [])
        location = json_dict.get('location')

        if not name:
            logger.error('Agent missing name.')
            return None
        if not db_refs:
            db_refs = {}
        agent = Agent(name, db_refs=db_refs)
        agent.mods = [ModCondition._from_json(mod) for mod in mods]
        agent.mutations = [MutCondition._from_json(mut) for mut in mutations]
        agent.bound_conditions = [BoundCondition._from_json(bc)
                                  for bc in bound_conditions]
        agent.location = location
        if activity:
            agent.activity = ActivityCondition._from_json(activity)
        return agent

    def __str__(self):
        attr_strs = []
        if self.mods:
            mod_str = 'mods: '
            mod_str += ', '.join(['%s' % m for m in self.mods])
            attr_strs.append(mod_str)
        if self.activity:
            if self.activity.is_active:
                attr_strs.append('%s' % self.activity.activity_type)
            else:
                attr_strs.append('%s: %s' % (self.activity.activity_type,
                                             self.activity.is_active))
        if self.mutations:
            mut_str = 'muts: '
            mut_str += ', '.join(['%s' % m for m in self.mutations])
            attr_strs.append(mut_str)
        if self.bound_conditions:
            attr_strs += ['bound: [%s, %s]' % (b.agent.name, b.is_bound)
                          for b in self.bound_conditions]
        if self.location:
            attr_strs += ['location: %s' % self.location]
        #if self.db_refs:
        #    attr_strs.append('db_refs: %s' % self.db_refs)
        attr_str = ', '.join(attr_strs)
        agent_name = self.name
        return '%s(%s)' % (agent_name, attr_str)


class BoundCondition(object):
    """Identify Agents bound (or not bound) to a given Agent in a given context.

    Parameters
    ----------
    agent : :py:class:`Agent`
        Instance of Agent.
    is_bound : bool
        Specifies whether the given Agent is bound or unbound in the current
        context. Default is True.

    Examples
    --------
    EGFR bound to EGF:

    >>> egf = Agent('EGF')
    >>> egfr = Agent('EGFR', bound_conditions=[BoundCondition(egf)])

    BRAF *not* bound to a 14-3-3 protein (YWHAB):

    >>> ywhab = Agent('YWHAB')
    >>> braf = Agent('BRAF', bound_conditions=[BoundCondition(ywhab, False)])
    """
    def __init__(self, agent, is_bound=True):
        self.agent = agent
        self.is_bound = is_bound

    def matches(self, other):
        return (self.matches_key() == other.matches_key())

    def matches_key(self):
        key = (self.agent.matches_key, self.is_bound)
        return str(key)

    def to_json(self):
        json_dict = _o({'agent': self.agent.to_json(),
                        'is_bound': self.is_bound})
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        agent_entry = json_dict.get('agent')
        if agent_entry is None:
            logger.error('BoundCondition missing agent.')
            return None
        agent = Agent._from_json(agent_entry)
        if agent is None:
            return None
        is_bound = json_dict.get('is_bound')
        if is_bound is None:
            logger.warning('BoundCondition missing is_bound, defaulting '
                           'to True.')
            is_bound = True
        bc = BoundCondition(agent, is_bound)
        return bc


class MutCondition(object):
    """Mutation state of an amino acid position of an Agent.

    Parameters
    ----------
    position : str
        Residue position of the mutation in the protein sequence.
    residue_from : str
        Wild-type (unmodified) amino acid residue at the given position.
    residue_to : str
        Amino acid at the position resulting from the mutation.

    Examples
    --------
    Represent EGFR with a L858R mutation:

    >>> egfr_mutant = Agent('EGFR', mutations=[MutCondition('858', 'L', 'R')])
    """
    def __init__(self, position, residue_from, residue_to=None):
        self.position = position
        self.residue_from = get_valid_residue(residue_from)
        self.residue_to = get_valid_residue(residue_to)

    def matches(self, other):
        return (self.matches_key() == other.matches_key())

    def matches_key(self):
        key = (str(self.position), str(self.residue_from),
               str(self.residue_to))
        return str(key)

    def equals(self, other):
        pos_match = (self.position == other.position)
        residue_from_match = (self.residue_from == other.residue_from)
        residue_to_match = (self.residue_to == other.residue_to)
        return (pos_match and residue_from_match and residue_to_match)

    def to_json(self):
        json_dict = _o({'position': self.position,
                        'residue_from': self.residue_from,
                        'residue_to': self.residue_to})
        return json_dict

    def to_hgvs(self):
        res_from = _aa_short_caps(self.residue_from)
        res_to = _aa_short_caps(self.residue_to)
        if res_to and res_from and self.position:
            hgvs_str = 'p.%s%s%s' % (res_from, self.position, res_to)
        elif res_to is None and res_from and self.position:
            hgvs_str = 'p.%s%s?' % (res_from, self.position)
        else:
            hgvs_str = 'p.?'
        return hgvs_str

    @classmethod
    def _from_json(cls, json_dict):
        position = json_dict.get('position')
        residue_from = json_dict.get('residue_from')
        residue_to = json_dict.get('residue_to')
        mc = cls(position, residue_from, residue_to)
        return mc

    def __str__(self):
        s = '(%s, %s, %s)' % (self.residue_from, self.position,
                              self.residue_to)
        return s

    def __repr__(self):
        return 'MutCondition' + str(self)

    def refinement_of(self, other):
        from_match = (self.residue_from == other.residue_from or
                      (self.residue_from is not None and other.residue_from is None))
        to_match = (self.residue_to == other.residue_to or
                    (self.residue_to is not None and other.residue_to is None))
        pos_match = (self.position == other.position or
                     (self.position is not None and other.position is None))
        return (from_match and to_match and pos_match)


class ModCondition(object):
    """Post-translational modification state at an amino acid position.

    Parameters
    ----------
    mod_type : str
        The type of post-translational modification, e.g., 'phosphorylation'.
        Valid modification types currently include: 'phosphorylation',
        'ubiquitination', 'sumoylation', 'hydroxylation', and 'acetylation'.
        If an invalid modification type is passed an InvalidModTypeError is
        raised.
    residue : str or None
        String indicating the modified amino acid, e.g., 'Y' or 'tyrosine'.
        If None, indicates that the residue at the modification site is
        unknown or unspecified.
    position : str or None
        String indicating the position of the modified amino acid, e.g., '202'.
        If None, indicates that the position is unknown or unspecified.
    is_modified : bool
        Specifies whether the modification is present or absent. Setting the
        flag specifies that the Agent with the ModCondition is unmodified
        at the site.

    Examples
    --------
    Doubly-phosphorylated MEK (MAP2K1):

    >>> phospho_mek = Agent('MAP2K1', mods=[
    ... ModCondition('phosphorylation', 'S', '202'),
    ... ModCondition('phosphorylation', 'S', '204')])

    ERK (MAPK1) unphosphorylated at tyrosine 187:

    >>> unphos_erk = Agent('MAPK1', mods=(
    ... ModCondition('phosphorylation', 'Y', '187', is_modified=False)))
    """
    def __init__(self, mod_type, residue=None, position=None,
                 is_modified=True):
        if mod_type not in modtype_conditions:
            logger.warning('Unknown modification type: %s' % mod_type)
        self.mod_type = mod_type
        self.residue = get_valid_residue(residue)
        if isinstance(position, int):
            self.position = str(position)
        else:
            self.position = position
        self.is_modified = is_modified

    def refinement_of(self, other, ontology):
        if self.is_modified != other.is_modified:
            return False
        type_match = (self.mod_type == other.mod_type or
                      ontology.isa('INDRA_MODS', self.mod_type,
                                   'INDRA_MODS', other.mod_type))
        residue_match = (self.residue == other.residue or
                         (self.residue is not None and other.residue is None))
        pos_match = (self.position == other.position or
                     (self.position is not None and other.position is None))
        return (type_match and residue_match and pos_match)

    def matches(self, other):
        return (self.matches_key() == other.matches_key())

    def matches_key(self):
        key = (str(self.mod_type), str(self.residue),
               str(self.position), str(self.is_modified))
        return str(key)

    def __str__(self):
        ms = '%s' % self.mod_type
        if self.residue is not None:
            ms += ', %s' % self.residue
        if self.position is not None:
            ms += ', %s' % self.position
        if not self.is_modified:
            ms += ', False'
        ms = '(' + ms + ')'
        return ms

    def __repr__(self):
        return str(self)

    def to_json(self):
        json_dict = _o({'mod_type': self.mod_type})
        if self.residue is not None:
            json_dict['residue'] = self.residue
        if self.position is not None:
            json_dict['position'] = self.position
        json_dict['is_modified'] = self.is_modified
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        mod_type = json_dict.get('mod_type')
        if not mod_type:
            logger.error('ModCondition missing mod_type.')
            return None
        if mod_type not in modtype_to_modclass.keys():
            logger.warning('Unknown modification type: %s' % mod_type)
        residue = json_dict.get('residue')
        position = json_dict.get('position')
        is_modified = json_dict.get('is_modified')
        if is_modified is None:
            logger.warning('ModCondition missing is_modified, defaulting '
                           'to True')
            is_modified = True
        mc = ModCondition(mod_type, residue, position, is_modified)
        return mc

    def equals(self, other):
        type_match = (self.mod_type == other.mod_type)
        residue_match = (self.residue == other.residue)
        pos_match = (self.position == other.position)
        is_mod_match = (self.is_modified == other.is_modified)
        return (type_match and residue_match and pos_match and is_mod_match)

    def __hash__(self):
        return hash(self.matches_key())


class ActivityCondition(object):
    """An active or inactive state of a protein.

    Examples
    --------
    Kinase-active MAP2K1:

    >>> mek_active = Agent('MAP2K1',
    ...                    activity=ActivityCondition('kinase', True))

    Transcriptionally inactive FOXO3:

    >>> foxo_inactive = Agent('FOXO3',
    ...                     activity=ActivityCondition('transcription', False))


    Parameters
    ----------
    activity_type : str
        The type of activity, e.g. 'kinase'. The basic, unspecified molecular
        activity is represented as 'activity'. Examples of other activity
        types are 'kinase', 'phosphatase', 'catalytic', 'transcription',
        etc.
    is_active : bool
        Specifies whether the given activity type is present or absent.
    """
    def __init__(self, activity_type, is_active):
        if activity_type not in activity_types:
            logger.warning('Invalid activity type: %s' % activity_type)
        self.activity_type = activity_type
        self.is_active = is_active

    def refinement_of(self, other, ontology):
        if self.is_active != other.is_active:
            return False
        if self.activity_type == other.activity_type:
            return True
        if ontology.isa('INDRA_ACTIVITIES', self.activity_type,
                        'INDRA_ACTIVITIES', other.activity_type):
            return True

    def equals(self, other):
        type_match = (self.activity_type == other.activity_type)
        is_act_match = (self.is_active == other.is_active)
        return (type_match and is_act_match)

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def matches_key(self):
        key = (str(self.activity_type), str(self.is_active))
        return str(key)

    def to_json(self):
        json_dict = _o({'activity_type': self.activity_type,
                        'is_active': self.is_active})
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        activity_type = json_dict.get('activity_type')
        is_active = json_dict.get('is_active')
        if not activity_type:
            logger.error('ActivityCondition missing activity_type, ' +
                         'defaulting to `activity`')
            activity_type = 'activity'
        if is_active is None:
            logger.warning('ActivityCondition missing is_active, ' +
                           'defaulting to True')
            is_active = True
        ac = ActivityCondition(activity_type, is_active)
        return ac

    def __str__(self):
        s = '%s' % self.activity_type
        if not self.is_active:
            s += ', False'
        s = '(' + s + ')'
        return s

    def __repr__(self):
        return str(self)


def _aa_short_caps(res):
    if res is None:
        return None
    res_info = amino_acids.get(res)
    if not res_info:
        return None
    return res_info['short_name'].capitalize()


def get_grounding(db_refs, ns_order=None):
    """Return a tuple of a preferred grounding namespace and ID.

    Parameters
    ----------
    db_refs : dict
        A dict of namespace to ID references associated with an agent.
    ns_order : list
        A list of namespaces which are in order of priority. The first
        matched namespace will be used as the grounding.

    Returns
    -------
    tuple
        A tuple whose first element is a grounding namespace (HGNC,
        CHEBI, etc.) and the second element is an identifier in the
        namespace. If no preferred grounding is available, a tuple of
        Nones is returned.
    """
    if ns_order is None:
        ns_order = default_ns_order
    for db_ns in ns_order:
        db_id = db_refs.get(db_ns)
        if not db_id:
            continue
        if isinstance(db_id, (list, tuple)):
            db_id = db_id[0]
        return db_ns, db_id
    return None, None
