import os
from collections import namedtuple
import textwrap

def read_amino_acids():
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

amino_acids, amino_acids_reverse = read_amino_acids()

BoundCondition = namedtuple('BoundCondition', ['agent', 'is_bound'])

class MutCondition(object):
    def __init__(self, position, residue_from, residue_to=None):
        self.position = position
        self.residue_from = get_valid_residue(residue_from)
        self.residue_to = get_valid_residue(residue_to)

    def matches(self, other):
        return (self.matches_key() == other.matches_key())

    def matches_key(self):
        key = (self.position, self.residue_from, self.residue_to)
        return str(key)

    def __eq__(self, other):
        pos_match = (self.position == other.position)
        residue_from_match = (self.residue_from == other.residue_from)
        residue_to_match = (self.residue_to == other.residue_to)
        return (pos_match and residue_from_match and residue_to_match)

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        s = '(%s, %s, %s)' % (self.residue_from, self.position,
                               self.residue_to)
        return s

class ModCondition(object):
    def __init__(self, mod_type, residue=None, position=None, is_modified=True):
        self.mod_type = mod_type
        self.residue = get_valid_residue(residue)
        if not isinstance(position, basestring):
            if position is None:
                self.position = None
            else:
                self.position = str(position)
        else:
            self.position = position
        self.is_modified = is_modified

    def refinement_of(self, other, mod_hierarchy):
        type_match = (self.mod_type == other.mod_type or \
            mod_hierarchy.isa(self.mod_type, other.mod_type))
        residue_match = (self.residue == other.residue or \
            (self.residue is not None and other.residue is None))
        pos_match = (self.position == other.position or \
            (self.position is not None and other.position is None))
        return (type_match and residue_match and pos_match)

    def matches(self, other):
        return (self.matches_key() == other.matches_key())

    def matches_key(self):
        key = (self.mod_type, self.residue, self.position, self.is_modified)
        return str(key)

    def __str__(self):
        ms = '%s' % self.mod_type
        if self.residue is not None:
            ms += ', %s' % self.residue
        if self.position is not None:
            ms += ', %s' % self.position
        if not self.is_modified:
            ms += ', FALSE'
        ms = '(' + ms + ')'
        return ms

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        type_match = (self.mod_type == other.mod_type)
        residue_match = (self.residue == other.residue)
        pos_match = (self.position == other.position)
        is_mod_match = (self.is_modified == other.is_modified)
        return (type_match and residue_match and pos_match and is_mod_match)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.matches_key())

class Agent(object):
    def __init__(self, name, mods=None, active=None,
                 bound_conditions=None, mutations=None, db_refs=None):
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

        if db_refs is None:
            self.db_refs = {}
        else:
            self.db_refs = db_refs

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def matches_key(self):
        # NOTE: Making a set of the mod matches_keys mat break if
        # you have an agent with two phosphorylations at serine
        # with unknown sites.
        key = (self.name,
               set([m.matches_key() for m in self.mods]),
               set([m.matches_key() for m in self.mutations]),
               self.active,
               len(self.bound_conditions),
               tuple((bc.agent.matches_key(), bc.is_bound)
                     for bc in sorted(self.bound_conditions,
                                      key=lambda x: x.agent.name)))
        return str(key)

    def entity_matches(self, other):
        return self.entity_matches_key() == other.entity_matches_key()

    def entity_matches_key(self):
        return self.name

    def refinement_of(self, other, entity_hierarchy, mod_hierarchy):
        # Make sure the Agent types match
        if type(self) != type(other):
            return False

        # ENTITIES
        # Check that the basic entity of the agent either matches or is related
        # to the entity of the other agent. If not, no match.
        if not (self.entity_matches(other) or \
                entity_hierarchy.isa(self.name, other.name)):
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
                    entity_hierarchy.isa(bc_self.agent.name,
                                         bc_other.agent.name)) and \
                    bc_self.is_bound == bc_other.is_bound:
                    bc_found = True
            # If we didn't find a match for this bound condition in self, then
            # no refinement
            if not bc_found:
                return False

        # MODIFICATIONS
        # Similar to the above, we check that self has all of the modifications
        # of other.
        # Make sure they have the same modifications
        for other_mod in other.mods:
            mod_found = False
            for self_mod in self.mods:
                if self_mod.refinement_of(other_mod, mod_hierarchy):
                    mod_found = True
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

        # Everything checks out
        return True

    def __eq__(self, other):
        matches = (self.name == other.name) and\
                  (self.mods == other.mods) and\
                  (self.mutations == other.mutations) and\
                  (self.bound_conditions == other.bound_conditions) and\
                  (self.active == other.active) and\
                  (self.db_refs == other.db_refs)
        return matches

    def __ne__(self, other):
        return not self == other

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
        #if self.db_refs:
        #    attr_strs.append('db_refs: %s' % self.db_refs)
        attr_str = ', '.join(attr_strs)
        return '%s(%s)' % (self.name, attr_str)

    def __repr__(self):
        return self.__str__()


class Evidence(object):
    """Container for evidence supporting a given statement.

    Attributes
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
    annotations : list
        List containing additional information on the context of the statement,
        e.g., species, cell line, tissue type, etc. The entries may vary
        depending on the source of the information.
    epistemics : string
        An identifier describing the epistemic certainty associated with the
        statement.
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
            self.annotations = []
        self.epistemics = epistemics

    def __eq__(self, other):
        matches = (self.source_api == other.source_api) and\
                  (self.source_id == other.source_id) and\
                  (self.pmid == other.pmid) and\
                  (self.text == other.text) and\
                  (self.annotations == other.annotations) and\
                  (self.epistemics == other.epistemics)
        return matches

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        ev_str = u'Evidence(%s, %s, %s, %s)' % \
                 (self.source_api, self.pmid, self.annotations,
                  self.text)
        return ev_str.encode('utf-8')

    def __repr__(self):
        return self.__str__()


class Statement(object):
    """The parent class of all statements.

    Attributes
    ----------
    evidence : list of Evidence objects.
        If a list of Evidence objects is passed to the constructor, the
        value is set to this list. If a bare Evidence object is passed,
        it is enclosed in a list. If no evidence is passed (the default),
        the value is set to an empty list.
    supports : list of Statements
        Statements that this Statement supports.
    supported_by : list of Statements
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

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def entities_match(self, other):
        return self.entities_match_key() == other.entities_match_key()

    def entities_match_key(self):
        key = (type(self), tuple(a.entity_matches_key() if a is not None
                                  else None for a in self.agent_list()))
        return str(key)

    def print_supports(self):
        print '%s supported_by:' % self.__str__()
        if self.supported_by:
            print '-->'
            for s in self.supported_by:
                s.print_supports()

    def __repr__(self):
        return self.__str__()


class Modification(Statement):
    """Generic statement representing the modification of a protein"""

    def __init__(self, enz, sub, residue=None, position=None, evidence=None):
        super(Modification, self).__init__(evidence)
        self.enz = enz
        self.sub = sub
        self.residue = get_valid_residue(residue)
        if position is not None:
            if not isinstance(position, basestring):
                position = str(position)
        self.position = position

    def matches_key(self):
        if self.enz is None:
            enz_key = None
        else:
            enz_key = self.enz.matches_key()
        key = (type(self), enz_key, self.sub.matches_key(),
               self.residue, self.position)
        return str(key)

    def agent_list(self):
        return [self.enz, self.sub]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Modification has two agents in agent_list.")
        self.enz = agent_list[0]
        self.sub = agent_list[1]

    def refinement_of(self, other, entity_hierarchy, mod_hierarchy):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if self.enz is None:
            enz_refinement = False
        elif self.enz is not None and other.enz is None:
            enz_refinement = True
        else:
            enz_refinement = self.enz.refinement_of(other.enz, entity_hierarchy,
                                                    mod_hierarchy)
        sub_refinement = self.sub.refinement_of(other.sub, entity_hierarchy,
                                                mod_hierarchy)
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

    def __str__(self):
        s = ("%s(%s, %s, %s, %s)" %
                  (type(self).__name__, self.enz, self.sub,
                   self.residue, self.position))
        return s


class SelfModification(Statement):
    """Generic statement representing the self modification of a protein"""

    def __init__(self, enz, residue=None, position=None, evidence=None):
        super(SelfModification, self).__init__(evidence)
        self.enz = enz
        self.residue = get_valid_residue(residue)
        if position is not None:
            if not isinstance(position, basestring):
                position = str(position)
        self.position = position

    def __str__(self):
        s = ("%s(%s, %s, %s)" %
             (type(self).__name__, self.enz.name,
              self.residue, self.position))
        return s

    def matches_key(self):
        key = (type(self), self.enz.matches_key(),
               self.residue, self.position)
        return str(key)

    def agent_list(self):
        return [self.enz]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("SelfModification has one agent.")
        self.enz = agent_list[0]

    def refinement_of(self, other, entity_hierarchy, mod_hierarchy):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.enz.refinement_of(other.enz, entity_hierarchy,
                                      mod_hierarchy):
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


class Phosphorylation(Modification):
    """Phosphorylation modification"""
    pass


class Autophosphorylation(SelfModification):
    """Autophosphorylation happens when a protein phosphorylates itself.

    A more precise name for this is cis-autophosphorylation.
    """
    pass


class Transphosphorylation(SelfModification):
    """Transphosphorylation assumes that a kinase is already bound to a
    substrate (usually of the same molecular species), and phosphorylates it in
    an intra-molecular fashion. The enz property of the statement must have
    exactly one bound_conditions entry, and we assume that enz phosphorylates
    this molecule. The bound_neg property is ignored here.  """
    pass


class Dephosphorylation(Modification):
    """Dephosphorylation modification"""
    pass


class Hydroxylation(Modification):
    """Hydroxylation modification"""
    pass


class Sumoylation(Modification):
    """Sumoylation modification"""
    pass


class Acetylation(Modification):
    """Acetylation modification"""
    pass


class Ubiquitination(Modification):
    """Ubiquitination modification"""
    pass


class ActivityActivity(Statement):
    """Statement representing the activation of a protein as a result of the
    activity of another protein."""

    def __init__(self, subj, subj_activity, relationship, obj,
                 obj_activity, evidence=None):
        super(ActivityActivity, self).__init__(evidence)
        self.subj = subj
        self.subj_activity = subj_activity
        self.obj = obj
        self.obj_activity = obj_activity
        self.relationship = relationship

    def matches_key(self):
        key = (type(self), self.subj.matches_key(), self.subj_activity,
                self.obj.matches_key(), self.obj_activity)
        return str(key)

    def agent_list(self):
        return [self.subj, self.obj]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("ActivityActivity has two agents.")
        self.subj = agent_list[0]
        self.obj = agent_list[1]

    def refinement_of(self, other, eh, mh):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        if self.subj.refinement_of(other.subj, eh, mh) and \
           self.obj.refinement_of(other.obj, eh, mh) and \
           self.subj_activity == other.subj_activity and \
           self.obj_activity == other.obj_activity and \
           self.relationship == other.relationship:
            return True
        else:
            return False

    def __str__(self):
        s = ("%s(%s, %s, %s, %s, %s)" %
             (type(self).__name__, self.subj, self.subj_activity,
              self.relationship, self.obj, self.obj_activity))
        return s


class RasGtpActivityActivity(ActivityActivity):
    pass


class ActiveForm(Statement):
    """Statement representing the activity of an Agent in a state described by
    one or more Agent conditions (including modification conditions, 
    bound-to conditions and mutation conditions)."""

    def __init__(self, agent, activity, is_active, evidence=None):
        super(ActiveForm, self).__init__(evidence)
        self.agent = agent
        self.activity = activity
        self.is_active = is_active

    def matches_key(self):
        key = (type(self), self.agent.matches_key(),
                self.activity, self.is_active)
        return str(key)

    def agent_list(self):
        return [self.agent]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("ActivityForm has one agent.")
        self.agent = agent_list[0]

    def refinement_of(self, other, entity_hierarchy, mod_hierarchy):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.agent.refinement_of(other.agent, entity_hierarchy,
                                          mod_hierarchy):
            return False

        # Make sure that the relationships and activities match
        # TODO: Develop an activity hierarchy? In which kinaseactivity is a
        # refinement of activity, for example.
        if self.activity == other.activity and \
           self.is_active == other.is_active:
               return True
        else:
            return False

    def __str__(self):
        s = ("ActiveForm(%s, %s, %s)" %
                (self.agent, self.activity, self.is_active))
        return s

class RasGef(Statement):
    """Statement representing the activation of a GTP-bound protein
    upon Gef activity."""

    def __init__(self, gef, gef_activity, ras, evidence=None):
        super(RasGef, self).__init__(evidence)
        self.gef = gef
        self.gef_activity = gef_activity
        self.ras = ras

    def matches_key(self):
        key = (type(self), self.gef.matches_key(), self.gef_activity,
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

    def refinement_of(self, other, eh, mh):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check the GEF
        if self.gef.refinement_of(other.gef, eh, mh) and \
           self.ras.refinement_of(other.ras, eh, mh) and \
           self.gef_activity == other.gef_activity:
            return True
        else:
            return False


class RasGap(Statement):
    """Statement representing the inactivation of a GTP-bound protein
    upon Gap activity."""

    def __init__(self, gap, gap_activity, ras, evidence=None):
        super(RasGap, self).__init__(evidence)
        self.gap = gap
        self.gap_activity = gap_activity
        self.ras = ras

    def matches_key(self):
        key = (type(self), self.gap.matches_key(), self.gap_activity,
                self.ras.matches_key())
        return str(key)

    def agent_list(self):
        return [self.gap, self.ras]

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("RasGap has two agents.")
        self.gap = agent_list[0]
        self.ras = agent_list[1]

    def refinement_of(self, other, eh, mh):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check the GAP
        if self.gap.refinement_of(other.gap, eh, mh) and \
           self.ras.refinement_of(other.ras, eh, mh) and \
           self.gap_activity == other.gap_activity:
            return True
        else:
            return False

    def __str__(self):
        s = ("RasGap(%s, %s, %s)" %
                (self.gap.name, self.gap_activity, self.ras.name))
        return s


class Complex(Statement):
    """Statement representing complex formation between a set of members"""

    def __init__(self, members, evidence=None):
        super(Complex, self).__init__(evidence)
        self.members = members

    def matches_key(self):
        key = (type(self), tuple(m.matches_key() for m in sorted(self.members,
                                                 key=lambda x: x.name)))
        return str(key)

    def entities_match_key(self):
        key = (type(self), tuple(a.entity_matches_key() if a is not None
                                  else None for a in sorted(self.members,
                                                key=lambda x: x.name)))
        return str(key)

    def agent_list(self):
        return self.members

    def set_agent_list(self, agent_list):
        self.members = agent_list

    def __str__(self):
        s = "Complex(%s)" % (', '.join([('%s' % m) for m in self.members]))
        return s

    def refinement_of(self, other, eh, mh):
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
                if self_agent.refinement_of(other_agent, eh, mh):
                    self_match_indices.add(self_agent_ix)
                    break
        if len(self_match_indices) != len(other.members):
            return False
        else:
            return True

def get_valid_residue(residue):
    if residue is not None and amino_acids.get(residue) is None:
        res = amino_acids_reverse.get(residue.lower())
        if res is None:
            raise InvalidResidueError(residue)
        else:
            return res
    return residue

class InvalidResidueError(ValueError):
    """Invalid residue (amino acid) name."""
    def __init__(self, name):
        ValueError.__init__(self, "Invalid residue name: '%s'" % name)
