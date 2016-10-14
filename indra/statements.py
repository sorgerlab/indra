"""
Statements representing mechanistic relationships between biological agents.

Statement classes follow an inheritance hierarchy, with all Statement types
inheriting from the parent class :py:class:`Statement`. At
the next level in the hierarchy are the following classes:

- :py:class:`Complex`
- :py:class:`Modification`
- :py:class:`SelfModification`
- :py:class:`RasGef`
- :py:class:`RasGap`
- :py:class:`Activation`

There are several types of Statements representing post-translational
modifications that further inherit from
:py:class:`Modification`:

- :py:class:`Phosphorylation`
- :py:class:`Dephosphorylation`
- :py:class:`Ubiquitination`
- :py:class:`Sumoylation`
- :py:class:`Hydroxylation`
- :py:class:`Acetylation`

There are additional subtypes of :py:class:`SelfModification`:

- :py:class:`Autophosphorylation`
- :py:class:`Transphosphorylation`

Statements involve one or more biological *Agents*, typically proteins,
represented by the class :py:class:`Agent`. Agents can have a specific
post-translational modification state (indicated by one or more instances of
:py:class:`ModCondition`), other bound Agents (:py:class:`BoundCondition`), or
mutations (:py:class:`MutCondition`). The *active* form of an agent (in terms
of its post-translational modifications or bound state) is indicated by an
instance of the class :py:class:`ActiveForm`.

Interactions between proteins are often described simply in terms of their
effect on a protein's "activity", e.g., "Active MEK activates ERK", or "DUSP6
inactives ERK".  These types of relationships are indicated by the statement
:py:class:`Activation`.

The evidence for a given Statement, which could include relevant citations,
database identifiers, and passages of text from the scientific literature, is
contained in one or more :py:class:`Evidence` objects associated with the
Statement.
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import os
import sys
import logging
import textwrap
import jsonpickle
from collections import namedtuple
import indra.databases.hgnc_client as hgc
import indra.databases.uniprot_client as upc

logger = logging.getLogger('indra_statements')

# Set the JSONpickle backend. We need to use the json module explicitly (rather
# than simplejson, which is the default if installed), because it returns
# unicode strings upon unpickling, which is what we want.
jsonpickle.set_preferred_backend('json')

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
    >>> egfr = Agent('EGFR', bound_conditions=(BoundCondition(egf)))

    BRAF *not* bound to a 14-3-3 protein (YWHAB):

    >>> ywhab = Agent('YWHAB')
    >>> braf = Agent('BRAF', bound_conditions=(BoundCondition(ywhab, False)))
    """
    def __init__(self, agent, is_bound=True):
        self.agent = agent
        self.is_bound = is_bound

@python_2_unicode_compatible
class MutCondition(object):
    """Mutation state of an amino acid position of an Agent.

    Parameters
    ----------
    position : string
        Residue position of the mutation in the protein sequence.
    residue_from : string
        Wild-type (unmodified) amino acid residue at the given position.
    residue_to : string
        Amino acid at the position resulting from the mutation.

    Examples
    --------
    Represent EGFR with a L858R mutation:

    >>> egfr_mutant = Agent('EGFR', mutations=(MutCondition('858', 'L', 'R')))
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

    def __str__(self):
        s = '(%s, %s, %s)' % (self.residue_from, self.position,
                               self.residue_to)
        return s

    def __repr__(self):
        return 'MutCondition' + str(self)


@python_2_unicode_compatible
class ModCondition(object):
    """Post-translational modification state at an amino acid position.

    Parameters
    ----------
    mod_type : string
        The type of post-translational modification, e.g., 'phosphorylation'.
        Valid modification types currently include: 'phosphorylation',
        'ubiquitination', 'sumoylation', 'hydroxylation', and 'acetylation'.
        If an invalid modification type is passed an InvalidModTypeError is
        raised.
    residue : string or None
        String indicating the modified amino acid, e.g., 'Y' or 'tyrosine'.
        If None, indicates that the residue at the modification site is
        unknown or unspecified.
    position : string or None
        String indicating the position of the modified amino acid, e.g., '202'.
        If None, indicates that the position is unknown or unspecified.
    is_modified : bool
        Specifies whether the modification is present or absent. Setting the
        flag specifies that the Agent with the ModCondition is unmodified
        at the site.

    Examples
    --------
    Doubly-phosphorylated MEK (MAP2K1):

    >>> phospho_mek = Agent('MAP2K1', mods=(
    ... ModCondition('phosphorylation', 'S', '202'),
    ... ModCondition('phosphorylation', 'S', '204')))

    ERK (MAPK1) unphosphorylated at tyrosine 187:

    >>> unphos_erk = Agent('MAPK1', mods=(
    ... ModCondition('phosphorylation', 'Y', '187', is_modified=False)))
    """
    def __init__(self, mod_type, residue=None, position=None, is_modified=True):
        self.mod_type = mod_type
        self.residue = get_valid_residue(residue)
        if isinstance(position, int):
            self.position = str(position)
        else:
            self.position = position
        self.is_modified = is_modified

    def refinement_of(self, other, mod_hierarchy):
        type_match = (self.mod_type == other.mod_type or \
            mod_hierarchy.isa('INDRA', self.mod_type, 'INDRA', other.mod_type))
        residue_match = (self.residue == other.residue or \
            (self.residue is not None and other.residue is None))
        pos_match = (self.position == other.position or \
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

    def equals(self, other):
        type_match = (self.mod_type == other.mod_type)
        residue_match = (self.residue == other.residue)
        pos_match = (self.position == other.position)
        is_mod_match = (self.is_modified == other.is_modified)
        return (type_match and residue_match and pos_match and is_mod_match)

    def __hash__(self):
        return hash(self.matches_key())

@python_2_unicode_compatible
class Agent(object):
    """A molecular entity, e.g., a protein.

    Parameters
    ----------
    name : string
        The name of the agent, preferably a canonicalized name such as an
        HGNC gene name.
    mods : list of :py:class:`ModCondition`
        Modification state of the agent.
    bound_conditions : list of :py:class:`BoundCondition`
        Other agents bound to the agent in this context.
    mutations : list of :py:class:`MutCondition`
        Amino acid mutations of the agent.
    location : str
        Cellular location of the agent. Must be a valid name (e.g. "nucleus")
        or identifier (e.g. "GO:0005634")for a GO cellular compartment.
    db_refs : dict
        Dictionary of database identifiers associated with this agent.
    """
    def __init__(self, name, mods=None, active=None,
                 bound_conditions=None, mutations=None,
                 location=None, db_refs=None):
        self.name = name

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

        self.active = active
        self.location = get_valid_location(location)

        if db_refs is None:
            self.db_refs = {}
        else:
            self.db_refs = db_refs

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def matches_key(self):
        # NOTE: Making a set of the mod matches_keys might break if
        # you have an agent with two phosphorylations at serine
        # with unknown sites.
        name_str = self.name
        key = (name_str,
               sorted([m.matches_key() for m in self.mods]),
               sorted([m.matches_key() for m in self.mutations]),
               self.active, self.location,
               len(self.bound_conditions),
               tuple((bc.agent.matches_key(), bc.is_bound)
                     for bc in sorted(self.bound_conditions,
                                      key=lambda x: x.agent.name)))
        return str(key)

    def entity_matches(self, other):
        return self.entity_matches_key() == other.entity_matches_key()

    def entity_matches_key(self):
        return self.name

    # Function to get the namespace to look in
    def get_grounding(self):
        be = self.db_refs.get('BE')
        if be:
            return ('BE', be)
        hgnc = self.db_refs.get('HGNC')
        if hgnc:
            if isinstance(hgnc, list):
                hgnc = hgnc[0]
            return ('HGNC', hgc.get_hgnc_name(str(hgnc)))
        up = self.db_refs.get('UP')
        if up:
            if isinstance(up, list):
                up = up[0]
            up_mnemonic = upc.get_mnemonic(up)
            if up_mnemonic and up_mnemonic.endswith('HUMAN'):
                gene_name = upc.get_gene_name(up, web_fallback=False)
                if gene_name:
                    return ('HGNC', gene_name)
            else:
                return ('UP', up)
        return (None, None)

    def isa(self, other, hierarchies):
        # Get the namespaces for the comparison
        (self_ns, self_id) = self.get_grounding()
        (other_ns, other_id) = other.get_grounding()
        # If one of the agents isn't grounded to a relevant namespace,
        # there can't be an isa relationship
        if not all((self_ns, self_id, other_ns, other_id)):
            return False
        # Check for isa relationship
        return hierarchies['entity'].isa(self_ns, self_id, other_ns, other_id)

    def refinement_of(self, other, hierarchies):
        # Make sure the Agent types match
        if type(self) != type(other):
            return False

        # ENTITIES
        # Check that the basic entity of the agent either matches or is related
        # to the entity of the other agent. If not, no match.

        # If the entities, match, then we can continue
        if not (self.entity_matches(other) or self.isa(other, hierarchies)):
            return False

        # BOUND CONDITIONS
        # Now check the bound conditions. For self to be a refinement of
        # other in terms of the bound conditions, it has to include all of the
        # bound conditions in the other agent, and add additional context.
        # TODO: For now, we do not check the bound conditions of the bound
        # conditions.
        # TODO: For now, we only look at exact agent matches, not at family
        # relationships among the bound conditions (this is to avoid the
        # confusion of relationships that might go in different directions
        # between the two statements).
        # FIXME: This matching procedure will get confused if the same
        # entity is included more than once in one of the sets--this will
        # be picked up as a match
        # Iterate over the bound conditions in the other agent, and make sure
        # they are all matched in self.
        for bc_other in other.bound_conditions:
            # Iterate over the bound conditions in self to find a match
            bc_found = False
            for bc_self in self.bound_conditions:
                if (bc_self.agent.entity_matches(bc_other.agent) or
                    bc_self.agent.isa(bc_other.agent, hierarchies)) and \
                    bc_self.is_bound == bc_other.is_bound:
                    bc_found = True
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
                if self_mod.refinement_of(other_mod,
                                          hierarchies['modification']):
                    # If this modification hasn't been used for matching yet
                    if not ix in matched_indices:
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
        # Make sure they have the same mutations
        for other_mut in other.mutations:
            mut_found = False
            for self_mut in self.mutations:
                if self_mut.matches(other_mut):
                    mut_found = True
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
            if not hierarchies['cellular_component'].partof(
                'INDRA', self.location, 'INDRA', other.location):
                return False

        # ACTIVITY
        if self.active is None:
            if other.active is not None:
                return False
        elif other.active is not None:
            if not hierarchies['activity'].isa('INDRA', self.active,
                                               'INDRA', other.active):
                return False

        # Everything checks out
        return True

    def equals(self, other):
        matches = (self.name == other.name) and\
                  (self.active == other.active) and\
                  (self.location == other.location) and\
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
                matches = matches and s.agent.equals(o.agent) and\
                                      s.is_bound == o.is_bound
        else:
            return False

        return matches

    def __str__(self):
        attr_strs = []
        if self.mods:
            mod_str = 'mods: '
            mod_str += ', '.join(['%s' % m for m in self.mods])
            attr_strs.append(mod_str)
        if self.active:
            attr_strs.append('active: %s' % self.active)
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

    def __repr__(self):
        return str(self)


@python_2_unicode_compatible
class Evidence(object):
    """Container for evidence supporting a given statement.

    Parameters
    ----------
    source_api : string or None
        String identifying the INDRA API used to capture the statement,
        e.g., 'trips', 'biopax', 'bel'.
    source_id : string or None
        For statements drawn from databases, ID of the database entity
        corresponding to the statement.
    pmid : string or None
        String indicating the Pubmed ID of the source of the statement.
    text : string
        Natural language text supporting the statement.
    annotations : dictionary
        Dictionary containing additional information on the
        context of the statement, e.g., species, cell line,
        tissue type, etc. The entries may vary depending on
        the source of the information.
    epistemics : dictionary
        A dictionary describing various forms of epistemic
        certainty associated with the statement.
    """
    def __init__(self, source_api=None, source_id=None, pmid=None, text=None,
                 annotations=None, epistemics=None):
        self.source_api = source_api
        self.source_id = source_id
        self.pmid = pmid
        self.text = text
        if annotations:
            self.annotations = annotations
        else:
            self.annotations = {}
        if epistemics:
            self.epistemics = epistemics
        else:
            self.epistemics = {}

    def equals(self, other):
        matches = (self.source_api == other.source_api) and\
                  (self.source_id == other.source_id) and\
                  (self.pmid == other.pmid) and\
                  (self.text == other.text) and\
                  (self.annotations == other.annotations) and\
                  (self.epistemics == other.epistemics)
        return matches

    def __str__(self):
        ev_str = 'Evidence(%s, %s, %s, %s)' % \
                 (self.source_api, self.pmid, self.annotations,
                  self.text)
        return ev_str

    def __repr__(self):
        if sys.version_info[0] >= 3:
            return str(self)
        else:
            return str(self).encode('utf-8')
        return str(self)


class Statement(object):
    """The parent class of all statements.

    Parameters
    ----------
    evidence : list of :py:class:`Evidence`
        If a list of Evidence objects is passed to the constructor, the
        value is set to this list. If a bare Evidence object is passed,
        it is enclosed in a list. If no evidence is passed (the default),
        the value is set to an empty list.
    supports : list of :py:class:`Statement`
        Statements that this Statement supports.
    supported_by : list of :py:class:`Statement`
        Statements supported by this statement.
    """
    def __init__(self, evidence=None, supports=None, supported_by=None):
        if evidence is None:
            self.evidence = []
        elif isinstance(evidence, Evidence):
            self.evidence = [evidence]
        elif isinstance(evidence, list):
            self.evidence = evidence
        else:
            raise ValueError('evidence must be an Evidence object, a list '
                             '(of Evidence objects), or None.')

        # Initialize supports/supported_by fields, which should be lists
        self.supports = supports if supports else []
        self.supported_by = supported_by if supported_by else []
        self.belief = 1

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def entities_match(self, other):
        return self.entities_match_key() == other.entities_match_key()

    def entities_match_key(self):
        key = (type(self), tuple(a.entity_matches_key() if a is not None
                                  else None for a in self.agent_list()))
        return str(key)

    def print_supports(self):
        print('%s supported_by:' % str(self))
        if self.supported_by:
            print('-->')
            for s in self.supported_by:
                s.print_supports()

    def __repr__(self):
        if sys.version_info[0] >= 3:
            return str(self)
        else:
            return str(self).encode('utf-8')

    def equals(self, other):
        if len(self.agent_list()) == len(other.agent_list()):
            for s, o in zip(self.agent_list(), other.agent_list()):
                if (s is None and o is not None) or\
                    (s is not None and o is None):
                    return False
                if s is not None and o is not None and not s.equals(o):
                    return False
        else:
            return False
        if len(self.evidence) == len(other.evidence):
            for s, o in zip(self.evidence, other.evidence):
                if not s.equals(o):
                    return False
        else:
            return False
        return True

    def to_json(self):
        """Return serialized Statement as a json string."""
        return jsonpickle.encode(self)

    @classmethod
    def from_json(cls, json_str):
        """Return Statement object by deserializing json string."""
        try:
            stmt = jsonpickle.decode(json_str)
            if isinstance(stmt, cls):
                return stmt
            else:
                logger.error('Could not construct Statement from json: ' + 
                             'Deserialized object is of type %s' % 
                             type(stmt).__name__)
                return None
        except ValueError as e:
            logger.error('Could not construct Statement from json: %s' % e)
            return None
        except IndexError as e:
            logger.error('Could not construct Statement from json: %s' % e)
            return None


@python_2_unicode_compatible
class Modification(Statement):
    """Generic statement representing the modification of a protein.

    Parameters
    ----------
    enz : :py:class`indra.statement.Agent`
        The enzyme involved in the modification.
    sub : :py:class:`indra.statement.Agent`
        The substrate of the modification.
    residue : string or None
        The amino acid residue being modified, or None if it is unknown or
        unspecified.
    position : string or None
        The position of the modified amino acid, or None if it is unknown or
        unspecified.
    evidence : list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """
    def __init__(self, enz, sub, residue=None, position=None, evidence=None):
        super(Modification, self).__init__(evidence)
        self.enz = enz
        self.sub = sub
        self.residue = get_valid_residue(residue)
        if isinstance(position, int):
            self.position = str(position)
        else:
            self.position = position

    def matches_key(self):
        if self.enz is None:
            enz_key = None
        else:
            enz_key = self.enz.matches_key()
        key = (type(self), enz_key, self.sub.matches_key(),
               str(self.residue), str(self.position))
        return str(key)

    def agent_list(self):
        return [self.enz, self.sub]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Modification has two agents in agent_list.")
        self.enz = agent_list[0]
        self.sub = agent_list[1]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if self.enz is None and other.enz is None:
            enz_refinement = True
        elif self.enz is None and other.enz is not None:
            enz_refinement = False
        elif self.enz is not None and other.enz is None:
            enz_refinement = True
        else:
            enz_refinement = self.enz.refinement_of(other.enz, hierarchies)
        sub_refinement = self.sub.refinement_of(other.sub, hierarchies)
        if not (enz_refinement and sub_refinement):
            return False
        # For this to be a refinement of the other, the modifications either
        # have to match or have this one be a subtype of the other; in
        # addition, the sites have to match, or this one has to have site
        # information and the other one not.
        residue_matches = (other.residue is None or\
                           (self.residue == other.residue))
        position_matches = (other.position is None or\
                            (self.position == other.position))
        return (residue_matches and position_matches)

    def equals(self, other):
        matches = super(Modification, self).equals(other)
        matches = matches and\
                  (self.residue == other.residue) and\
                  (self.position == other.position)
        return matches

    def __str__(self):
        res_str = (', %s' % self.residue) if self.residue is not None else ''
        pos_str = (', %s' % self.position) if self.position is not None else ''
        s = ("%s(%s, %s%s%s)" %
                  (type(self).__name__, self.enz, self.sub,
                   res_str, pos_str))
        return s


@python_2_unicode_compatible
class SelfModification(Statement):
    """Generic statement representing the self-modification of a protein.

    Parameters
    ----------
    enz : :py:class`indra.statement.Agent`
        The enzyme involved in the modification, which is also the substrate.
    residue : string or None
        The amino acid residue being modified, or None if it is unknown or
        unspecified.
    position : string or None
        The position of the modified amino acid, or None if it is unknown or
        unspecified.
    evidence : list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """
    def __init__(self, enz, residue=None, position=None, evidence=None):
        super(SelfModification, self).__init__(evidence)
        self.enz = enz
        self.residue = get_valid_residue(residue)
        if isinstance(position, int):
            self.position = str(position)
        else:
            self.position = position

    def __str__(self):
        res_str = (', %s' % self.residue) if self.residue is not None else ''
        pos_str = (', %s' % self.position) if self.position is not None else ''
        s = ("%s(%s%s%s)" %
                  (type(self).__name__, self.enz,
                   res_str, pos_str))
        return s

    def matches_key(self):
        key = (type(self), self.enz.matches_key(),
               str(self.residue), str(self.position))
        return str(key)

    def agent_list(self):
        return [self.enz]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("SelfModification has one agent.")
        self.enz = agent_list[0]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.enz.refinement_of(other.enz, hierarchies):
            return False
        # For this to be a refinement of the other, the modifications either
        # have to match or have this one be a subtype of the other; in
        # addition, the sites have to match, or this one has to have site
        # information and the other one not.
        residue_matches = (other.residue is None or\
                           (self.residue == other.residue))
        position_matches = (other.position is None or\
                            (self.position == other.position))
        return (residue_matches and position_matches)

    def equals(self, other):
        matches = super(SelfModification, self).equals(other)
        matches = matches and\
                  (self.residue == other.residue) and\
                  (self.position == other.position)
        return matches


class Phosphorylation(Modification):
    """Phosphorylation modification.

    Examples
    --------
    MEK (MAP2K1) phosphorylates ERK (MAPK1) at threonine 185:

    >>> mek = Agent('MAP2K1')
    >>> erk = Agent('MAPK1')
    >>> phos = Phosphorylation(mek, erk, 'T', '185')
    """
    pass


class Autophosphorylation(SelfModification):
    """Intramolecular autophosphorylation, i.e., in *cis*.

    Examples
    --------
    p38 bound to TAB1 cis-autophosphorylates itself (see :pmid:`19155529`).

    >>> tab1 = Agent('TAB1')
    >>> p38_tab1 = Agent('P38', bound_conditions=(BoundCondition(tab1)))
    >>> autophos = Autophosphorylation(p38_tab1)
    """
    pass


class Transphosphorylation(SelfModification):
    """Autophosphorylation in *trans.*

    Transphosphorylation assumes that a kinase is already bound to a substrate
    (usually of the same molecular species), and phosphorylates it in an
    intra-molecular fashion. The enz property of the statement must have
    exactly one bound_conditions entry, and we assume that enz phosphorylates
    this molecule. The bound_neg property is ignored here.
    """
    pass


class Dephosphorylation(Modification):
    """Dephosphorylation modification.

    Examples
    --------
    DUSP6 dephosphorylates ERK (MAPK1) at T185:

    >>> dusp6 = Agent('DUSP6')
    >>> erk = Agent('MAPK1')
    >>> dephos = Dephosphorylation(dusp6, erk, 'T', '185')
    """
    pass


class Hydroxylation(Modification):
    """Hydroxylation modification."""
    pass


class Dehydroxylation(Modification):
    """Dehydroxylation modification."""
    pass


class Sumoylation(Modification):
    """Sumoylation modification."""
    pass


class Desumoylation(Modification):
    """Desumoylation modification."""
    pass


class Acetylation(Modification):
    """Acetylation modification."""
    pass


class Deacetylation(Modification):
    """Deacetylation modification."""
    pass


class Glycosylation(Modification):
    """Glycosylation modification."""
    pass


class Deglycosylation(Modification):
    """Deglycosylation modification."""
    pass


class Ubiquitination(Modification):
    """Ubiquitination modification."""
    pass


class Deubiquitination(Modification):
    """Deubiquitination modification."""
    pass


class Farnesylation(Modification):
    """Farnesylation modification."""
    pass

@python_2_unicode_compatible
class Activation(Statement):
    """Indicates that the activity of a protein affects the activity of another.

    This statement is intended to be used for physical interactions where the
    mechanism of activation is not explicitly specified, which is often the
    case for descriptions of mechanisms extracted from the literature. Both
    activating and inactivating interactions can be represented.

    Parameters
    ----------
    subj : :py:class:`Agent`
        The agent responsible for the change in activity, i.e., the "upstream"
        node.
    subj_activity : string
        The type of biochemical activity responsible for the effect, e.g.,
        the subject's "kinase" activity.
    obj : :py:class:`Agent`
        The agent whose activity is influenced by the subject, i.e., the
        "downstream" node.
    obj_activity : string
        The activity of the obj Agent that is affected, e.g., its "kinase"
        activity.
    is_activation : bool
        Indicates the type of interaction: True for activation and
        False for inactivation/inhibition
    evidence : list of :py:class:`Evidence`
        Evidence objects in support of the modification.

    Examples
    --------

    The kinase activity of MEK (MAP2K1) activates the kinase activity of ERK
    (MAPK1):

    >>> mek = Agent('MAP2K1')
    >>> erk = Agent('MAPK1')
    >>> act = Activation(mek, 'kinase', erk, 'kinase', True)
    """
    def __init__(self, subj, subj_activity, obj, obj_activity, is_activation,
                 evidence=None):
        super(Activation, self).__init__(evidence)
        self.subj = subj
        self.subj_activity = subj_activity
        self.obj = obj
        self.obj_activity = obj_activity
        self.is_activation = is_activation

    def matches_key(self):
        key = (type(self), self.subj.matches_key(), str(self.subj_activity),
                self.obj.matches_key(), str(self.obj_activity))
        return str(key)

    def agent_list(self):
        return [self.subj, self.obj]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Activation has two agents.")
        self.subj = agent_list[0]
        self.obj = agent_list[1]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        if self.subj.refinement_of(other.subj, hierarchies) and \
            self.obj.refinement_of(other.obj, hierarchies):
            if self.is_activation != other.is_activation:
                return False
            subj_act_match = (self.subj_activity == other.subj_activity) or \
                hierarchies['activity'].isa('INDRA', self.subj_activity,
                                            'INDRA', other.subj_activity)
            obj_act_match = (self.obj_activity == other.obj_activity) or \
                hierarchies['activity'].isa('INDRA', self.obj_activity,
                                            'INDRA', other.obj_activity)
            if subj_act_match and obj_act_match:
                return True
            else:
                return False
        else:
            return False

    def __str__(self):
        s = ("%s(%s, %s, %s, %s, %s)" %
             (type(self).__name__, self.subj, self.subj_activity,
              self.obj, self.obj_activity, self.is_activation))
        return s

    def equals(self, other):
        matches = super(Activation, self).equals(other)
        matches = matches and\
                  (self.subj_activity == other.subj_activity) and\
                  (self.obj_activity == other.obj_activity) and\
                  (self.is_activation == other.is_activation)
        return matches


class RasGtpActivation(Activation):
    pass


@python_2_unicode_compatible
class ActiveForm(Statement):
    """Specifies conditions causing an Agent to be active or inactive.

    Types of conditions influencing a specific type of biochemical activity can
    include modifications, bound Agents, and mutations.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The Agent in a particular active or inactive state. The sets
        of ModConditions, BoundConditions, and MutConditions on the given
        Agent instance indicate the relevant conditions.
    activity : string
        The type of activity influenced by the given set of conditions,
        e.g., "kinase".
    is_active : bool
        Whether the conditions are activating (True) or inactivating (False).
    """
    def __init__(self, agent, activity, is_active, evidence=None):
        super(ActiveForm, self).__init__(evidence)
        self.agent = agent
        self.activity = activity
        self.is_active = is_active

    def matches_key(self):
        key = (type(self), self.agent.matches_key(),
                str(self.activity), str(self.is_active))
        return str(key)

    def agent_list(self):
        return [self.agent]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("ActiveForm has one agent.")
        self.agent = agent_list[0]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.agent.refinement_of(other.agent, hierarchies):
            return False

        # Make sure that the relationships and activities match
        if (self.is_active == other.is_active) and \
            (self.activity == other.activity or \
            hierarchies['activity'].isa('INDRA', self.activity,
                                        'INDRA', other.activity)):
               return True
        else:
            return False

    def __str__(self):
        s = ("ActiveForm(%s, %s, %s)" %
                (self.agent, self.activity, self.is_active))
        return s

    def equals(self, other):
        matches = super(ActiveForm, self).equals(other)
        matches = matches and\
                  (self.activity == other.activity) and\
                  (self.is_active == other.is_active)
        return matches


@python_2_unicode_compatible
class HasActivity(Statement):
    """States that an Agent has or doesn't have a given activity type.

    With this Statement, one cane express that a given protein is a kinase, or,
    for instance, that it is a transcription factor. It is also possible to
    construct negative statements with which one epxresses, for instance,
    that a given protein is not a kinase.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The Agent that that statement is about. Note that the detailed state
        of the Agent is not relevant for this type of statement.
    activity : string
        The type of activity, e.g., "kinase".
    has_activity : bool
        Whether the given Agent has the given activity (True) or
        not (False).
    """
    def __init__(self, agent, activity, has_activity, evidence=None):
        super(HasActivity, self).__init__(evidence)
        self.agent = agent
        self.activity = activity
        self.has_activity = has_activity

    def matches_key(self):
        key = (type(self), self.agent.matches_key(),
                str(self.activity), str(self.has_activity))
        return str(key)

    def agent_list(self):
        return [self.agent]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("HasActivity has one agent.")
        self.agent = agent_list[0]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.agent.refinement_of(other.agent, hierarchies):
            return False

        # Make sure that the relationships and activities match
        if (self.has_activity == other.has_activity) and \
            (self.activity == other.activity or \
            hierarchies['activity'].isa(self.activity, other.activity)):
               return True
        else:
            return False

    def __str__(self):
        s = ("HasActivity(%s, %s, %s)" %
                (self.agent, self.activity, self.has_activity))
        return s

    def equals(self, other):
        matches = super(HasActivity, self).equals(other)
        matches = matches and\
                  (self.activity == other.activity) and\
                  (self.has_activity == other.has_activity)
        return matches


@python_2_unicode_compatible
class RasGef(Statement):
    """Exchange of GTP for GDP on a Ras-family protein mediated by a GEF.

    Represents the generic process by which a guanosine exchange factor (GEF)
    catalyzes nucleotide exchange on a particular Ras superfamily protein.

    Parameters
    ----------
    gef : :py:class:`Agent`
        The guanosine exchange factor.
    gef_activity : string
        The biochemical activity of the GEF responsible for exchange.
    ras : :py:class:`Agent`
        The Ras superfamily protein.

    Examples
    --------
    SOS1 catalyzes nucleotide exchange on KRAS:

    >>> sos = Agent('SOS1')
    >>> kras = Agent('KRAS')
    >>> rasgef = RasGef(sos, 'gef', kras)
    """
    def __init__(self, gef, gef_activity, ras, evidence=None):
        super(RasGef, self).__init__(evidence)
        self.gef = gef
        self.gef_activity = gef_activity
        self.ras = ras

    def matches_key(self):
        key = (type(self), self.gef.matches_key(), str(self.gef_activity),
                self.ras.matches_key())
        return str(key)

    def agent_list(self):
        return [self.gef, self.ras]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("RasGef has two agents.")
        self.gef = agent_list[0]
        self.ras = agent_list[1]

    def __str__(self):
        s = ("RasGef(%s, %s, %s)" %
                (self.gef.name, self.gef_activity, self.ras.name))
        return s

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check the GEF
        if self.gef.refinement_of(other.gef, hierarchies) and \
           self.ras.refinement_of(other.ras, hierarchies) and \
           self.gef_activity == other.gef_activity:
            return True
        else:
            return False

    def equals(self, other):
        matches = super(RasGef, self).equals(other)
        matches = matches and (self.gef_activity == other.gef_activity)
        return matches


@python_2_unicode_compatible
class RasGap(Statement):
    """Acceleration of a Ras protein's GTP hydrolysis rate by a GAP.

    Represents the generic process by which a GTPase activating protein (GAP)
    catalyzes GTP hydrolysis by a particular Ras superfamily protein.

    Parameters
    ----------
    gap : :py:class:`Agent`
        The GTPase activating protein.
    gap_activity : string
        The biochemical activity of the GAP responsible for hydrolysis.
    ras : :py:class:`Agent`
        The Ras superfamily protein.

    Examples
    --------
    RASA1 catalyzes GTP hydrolysis on KRAS:

    >>> rasa1 = Agent('RASA1')
    >>> kras = Agent('KRAS')
    >>> rasgap = RasGap(rasa1, 'gap', kras)
    """
    def __init__(self, gap, gap_activity, ras, evidence=None):
        super(RasGap, self).__init__(evidence)
        self.gap = gap
        self.gap_activity = gap_activity
        self.ras = ras

    def matches_key(self):
        key = (type(self), self.gap.matches_key(), str(self.gap_activity),
                self.ras.matches_key())
        return str(key)

    def agent_list(self):
        return [self.gap, self.ras]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("RasGap has two agents.")
        self.gap = agent_list[0]
        self.ras = agent_list[1]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check the GAP
        if self.gap.refinement_of(other.gap, hierarchies) and \
           self.ras.refinement_of(other.ras, hierarchies) and \
           self.gap_activity == other.gap_activity:
            return True
        else:
            return False

    def __str__(self):
        s = ("RasGap(%s, %s, %s)" %
                (self.gap.name, self.gap_activity, self.ras.name))
        return s

    def equals(self, other):
        matches = super(RasGap, self).equals(other)
        matches = matches and (self.gap_activity == other.gap_activity)
        return matches


@python_2_unicode_compatible
class Complex(Statement):
    """A set of proteins observed to be in a complex.

    Parameters
    ----------
    members : list of :py:class:`Agent`
        The set of proteins in the complex.

    Examples
    --------
    BRAF is observed to be in a complex with RAF1:

    >>> braf = Agent('BRAF')
    >>> raf1 = Agent('RAF1')
    >>> cplx = Complex([braf, raf1])
    """
    def __init__(self, members, evidence=None):
        super(Complex, self).__init__(evidence)
        self.members = members

    def matches_key(self):
        key = (type(self), tuple(m.matches_key() for m in sorted(self.members,
                                                 key=lambda x: x.matches_key())))
        return str(key)

    def entities_match_key(self):
        key = (type(self), tuple(a.entity_matches_key() if a is not None
                                  else None for a in sorted(self.members,
                                                key=lambda x: x.matches_key())))
        return str(key)

    def agent_list(self):
        return self.members

    def set_agent_list(self, agent_list):
        self.members = agent_list

    def __str__(self):
        s = "Complex(%s)" % (', '.join([('%s' % m) for m in self.members]))
        return s

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Make sure the length of the members list is the same. Note that this
        # treats Complex([A, B, C]) as distinct from Complex([A, B]), rather
        # than as a refinement.
        if len(self.members) != len(other.members):
            return False
        # Check that every member in other is refined in self, but only once!
        self_match_indices = set([])
        for other_agent in other.members:
            for self_agent_ix, self_agent in enumerate(self.members):
                if self_agent_ix in self_match_indices:
                    continue
                if self_agent.refinement_of(other_agent, hierarchies):
                    self_match_indices.add(self_agent_ix)
                    break
        if len(self_match_indices) != len(other.members):
            return False
        else:
            return True

    def equals(self, other):
        matches = super(Complex, self).equals(other)
        return matches


@python_2_unicode_compatible
class Translocation(Statement):
    """The translocation of a molecular agent from one location to another.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The agent which translocates.
    from_location : Optional[str]
        The location from which the agent translocates. This must
        be a valid GO cellular component name (e.g. "cytoplasm")
        or ID (e.g. "GO:0005737").
    to_location : Optional[str]
        The location to which the agent translocates. This must
        be a valid GO cellular component name or ID.
    """
    def __init__(self, agent, from_location=None, to_location=None,
                 evidence=None):
        super(Translocation, self).__init__(evidence)
        self.agent = agent
        self.from_location = get_valid_location(from_location)
        self.to_location = get_valid_location(to_location)

    def agent_list(self):
        return [self.agent]

    def set_agent_list(self, agent_list):
        if(len(agent_list) != 1):
            raise ValueError("Translocation has 1 agent")
        self.agent = agent_list[0]

    def __str__(self):
        s = ("Translocation(%s, %s, %s)" %
                (self.agent, self.from_location, self.to_location))
        return s

    def refinement_of(self, other, hierarchies=None):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check several conditions for refinement
        ch = hierarchies['cellular_component']
        ref1 = self.agent.refinement_of(other.agent, hierarchies)
        ref2 = (other.from_location is None or
                self.from_location == other.from_location or
                ch.partof('INDRA', self.from_location,
                          'INDRA', other.from_location))
        ref3 = (other.to_location is None or
                self.to_location == other.to_location or
                ch.partof('INDRA', self.to_location,
                          'INDRA', other.to_location))
        return (ref1 and ref2 and ref3)

    def equals(self, other):
        matches = super(Translocation, self).equals(other)
        matches = matches and (self.from_location == other.from_location)
        matches = matches and (self.to_location == other.to_location)
        return matches

    def matches_key(self):
        key = (type(self), self.agent.matches_key(), str(self.from_location),
                str(self.to_location))
        return str(key)


def get_valid_residue(residue):
    """Check if the given string represents a valid amino acid residue."""
    if residue is not None and amino_acids.get(residue) is None:
        res = amino_acids_reverse.get(residue.lower())
        if res is None:
            raise InvalidResidueError(residue)
        else:
            return res
    return residue


def get_valid_location(location):
    """Check if the given location represents a valid cellular component."""
    # If we're given None, return None
    if location is not None and cellular_components.get(location) is None:
        loc = cellular_components_reverse.get(location)
        if loc is None:
            raise InvalidLocationError(location)
        else:
            return loc
    return location


def _read_cellular_components():
    """Read cellular components from a resource file."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    cc_file = this_dir + '/resources/cellular_components.tsv'
    cellular_components = {}
    cellular_components_reverse = {}
    with open(cc_file, 'rt') as fh:
        lines = fh.readlines()
    for lin in lines[1:]:
        terms = lin.strip().split('\t')
        cellular_components[terms[1]] = terms[0]
        cellular_components_reverse[terms[0]] = terms[1]
    return cellular_components, cellular_components_reverse


cellular_components, cellular_components_reverse = _read_cellular_components()


def _read_amino_acids():
    """Read the amino acid information from a resource file."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    aa_file = this_dir + '/resources/amino_acids.tsv'
    amino_acids = {}
    amino_acids_reverse = {}
    with open(aa_file, 'rt') as fh:
        lines = fh.readlines()
    for lin in lines[1:]:
        terms = lin.strip().split('\t')
        key = terms[2]
        val = {'full_name': terms[0],
               'short_name': terms[1],
               'indra_name': terms[3]}
        amino_acids[key] = val
        for v in val.values():
            amino_acids_reverse[v] = key
    return amino_acids, amino_acids_reverse

amino_acids, amino_acids_reverse = _read_amino_acids()


class InvalidResidueError(ValueError):
    """Invalid residue (amino acid) name."""
    def __init__(self, name):
        ValueError.__init__(self, "Invalid residue name: '%s'" % name)


class InvalidLocationError(ValueError):
    """Invalid cellular component name."""
    def __init__(self, name):
        ValueError.__init__(self, "Invalid location name: '%s'" % name)


