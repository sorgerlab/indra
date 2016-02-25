from pysb import *
from pysb import ReactionPattern, ComplexPattern, ComponentDuplicateNameError
from collections import namedtuple
import textwrap

BoundCondition = namedtuple('BoundCondition', ['agent', 'is_bound'])

class Agent(object):

    def __init__(self, name, mods=None, mod_sites=None, active=None,
                 bound_conditions=None, db_refs=None):
        self.name = name

        if mods is None:
            self.mods = []
        else:
            self.mods = mods

        if mod_sites is None:
            self.mod_sites = []
        else:
            self.mod_sites = mod_sites

        if bound_conditions is None:
            self.bound_conditions = []
        else:
            self.bound_conditions = bound_conditions

        self.active = active

        if db_refs is None:
            self.db_refs = {}
        else:
            self.db_refs = db_refs

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def matches_key(self):
        key = (self.name,
               set(self.mods),
               set(self.mod_sites),
               self.active,
               len(self.bound_conditions),
               tuple((bc.agent.matches_key(), bc.is_bound)
                     for bc in sorted(self.bound_conditions,
                                      key=lambda x: x.agent.name)))
        return key

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
        assert len(self.mods) == len(self.mod_sites), \
               "Mods and mod_sites must match."
        assert len(other.mods) == len(other.mod_sites), \
               "Mods and mod_sites must match."
        # Make sure they have the same modifications
        for other_mod_ix in range(len(other.mods)):
            mod_found = False
            other_mod = other.mods[other_mod_ix]
            other_mod_site = other.mod_sites[other_mod_ix]
            for self_mod_ix in range(len(self.mods)):
                self_mod = self.mods[self_mod_ix]
                self_mod_site = self.mod_sites[self_mod_ix]
                # Or has an isa relationship...
                if (self_mod == other_mod or \
                        mod_hierarchy.isa(self_mod, other_mod)) and \
                   (self_mod_site == other_mod_site or \
                        (self_mod_site is not None and other_mod_site is None)):
                    mod_found = True
            # If we didn't find an exact match for this mod in other, then
            # no refinement
            if not mod_found:
                return False

        # Everything checks out
        return True

    def __repr__(self):
        attr_strs = []
        if self.mods:
            attr_strs.append('mods: %s' % self.mods)
        if self.mod_sites:
            attr_strs.append('mod_sites: %s' % self.mod_sites)
        if self.active:
            attr_strs.append('active: %s' % self.active)
        if self.bound_conditions:
            attr_strs += ['bound: [%s, %s]' % (b.agent.name, b.is_bound)
                          for b in self.bound_conditions]
        if self.db_refs:
            attr_strs.append('db_refs: %s' % self.db_refs)
        attr_str = ', '.join(attr_strs)
        return '%s(%s)' % (self.name, attr_str)


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
        return (type(self), tuple(a.entity_matches_key()
                                  for a in self.agent_list()))

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

    def __init__(self, enz, sub, mod, mod_pos, evidence=None):
        super(Modification, self).__init__(evidence)
        self.enz = enz
        self.sub = sub
        self.mod = mod
        self.mod_pos = mod_pos

    def matches_key(self):
        return (type(self), self.enz.matches_key(), self.sub.matches_key(),
                self.mod, self.mod_pos)

    def agent_list(self):
        return [self.enz, self.sub]

    def refinement_of(self, other, entity_hierarchy, mod_hierarchy):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not (self.enz.refinement_of(other.enz, entity_hierarchy,
                                  mod_hierarchy) and \
                self.sub.refinement_of(other.sub, entity_hierarchy,
                                  mod_hierarchy)):
            return False
        # For this to be a refinement of the other, the modifications either
        # have to match or have this one be a subtype of the other; in
        # addition, the sites have to match, or this one has to have site
        # information and the other one not.
        if (self.mod == other.mod or \
                mod_hierarchy.isa(self.mod, other.mod)) and \
           (self.mod_pos == other.mod_pos or \
                (self.mod_pos is not None and other.mod_pos is None)):
            return True
        else:
            return False

    def __str__(self):
        s = ("%s(%s, %s, %s, %s)" %
                  (type(self).__name__, self.enz.name, self.sub.name, self.mod,
                   self.mod_pos))
        return s


class SelfModification(Statement):
    """Generic statement representing the self modification of a protein"""

    def __init__(self, enz, mod, mod_pos, evidence=None):
        super(SelfModification, self).__init__(evidence)
        self.enz = enz
        self.mod = mod
        self.mod_pos = mod_pos

    def __str__(self):
        s = ("%s(%s, %s, %s)" %
             (type(self).__name__, self.enz.name, self.mod, self.mod_pos))
        return s

    def matches_key(self):
        return (type(self), self.enz.matches_key(), self.mod, self.mod_pos)

    def agent_list(self):
        return [self.enz]

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
        if (self.mod == other.mod or \
                mod_hierarchy.isa(self.mod, other.mod)) and \
           (self.mod_pos == other.mod_pos or \
                (self.mod_pos is not None and other.mod_pos is None)):
            return True
        else:
            return False


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
        return (type(self), self.subj.matches_key(), self.subj_activity,
                self.obj.matches_key(), self.obj_activity)

    def agent_list(self):
        return [self.subj, self.obj]

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
             (type(self).__name__, self.subj.name, self.subj_activity,
              self.relationship, self.obj.name, self.obj_activity))
        return s


class RasGtpActivityActivity(ActivityActivity):
    pass


class ActivityModification(Statement):
    """Statement representing the activation of a protein as a result
    of a residue modification"""

    def __init__(self, monomer, mod, mod_pos, relationship, activity,
                 evidence=None):
        super(ActivityModification, self).__init__(evidence)
        self.monomer = monomer
        # This means that one of mod and mod_pos are a list, but not both,
        # which should raise a ValueError
        if isinstance(mod, list) and len(mod) > 1 and mod_pos is None:
            raise ValueError('If more than one modification is specified then '
                             'mod_pos must also be specified and cannot be '
                             'None.')
        elif isinstance(mod, list) != isinstance(mod_pos, list):
            raise ValueError('If mod or mod_pos are provided as lists they '
                             'must be matched.')
        elif isinstance(mod, list) and isinstance(mod_pos, list) and \
                len(mod) != len(mod_pos):
            raise ValueError('If mod and mod_pos are lists, then they must be '
                             'the same length.')
        elif isinstance(mod, list) and isinstance(mod_pos, list) and \
                len(mod) == len(mod_pos):
            # Set the fields to be the lists, but make sure to stringify
            # any entries in the mod_pos list that might be ints
            self.mod = mod
            self.mod_pos = [str(mp) for mp in mod_pos]
        elif isinstance(mod, basestring) and \
                (isinstance(mod_pos, int) or isinstance(mod_pos, basestring)):
            self.mod = [mod]
            self.mod_pos = [str(mod_pos)]
        elif isinstance(mod, basestring) and mod_pos is None:
            self.mod = [mod]
            self.mod_pos = [None]
        else:
            raise ValueError('Invalid values for mod and/or mod_pos')

        self.relationship = relationship
        self.activity = activity

    def matches_key(self):
        return (type(self), self.monomer.matches_key(), self.mod, self.mod_pos,
                self.relationship, self.activity)

    def agent_list(self):
        return [self.monomer]

    def refinement_of(self, other, entity_hierarchy, mod_hierarchy):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.monomer.refinement_of(other.monomer, entity_hierarchy,
                                          mod_hierarchy):
            return False

        # Mod and mod_pos should always be lists of the same length; mod_pos
        # can also be None
        assert isinstance(self.mod, list)
        assert isinstance(self.mod_pos, list) or mod_pos is None
        assert isinstance(other.mod, list)
        assert isinstance(other.mod_pos, list) or mod_pos is None
        # Make sure that every instance of a modification in other is also
        # found (or refined) in self. To facilitate comparisons, we first zip
        # the two lists together in a list of tuples that can be sorted to
        # canonicalize the ordering.

        # Make sure they have the same modifications
        for other_mod_ix in range(len(other.mod)):
            mod_found = False
            other_mod = other.mod[other_mod_ix]
            other_mod_pos = other.mod_pos[other_mod_ix]
            for self_mod_ix in range(len(self.mod)):
                self_mod = self.mod[self_mod_ix]
                self_mod_pos = self.mod_pos[self_mod_ix]
                # Or has an isa relationship...
                if (self_mod == other_mod or \
                        mod_hierarchy.isa(self_mod, other_mod)) and \
                   (self_mod_pos == other_mod_pos or \
                        (self_mod_pos is not None and other_mod_pos is None)):
                    mod_found = True
            # If we didn't find an exact match for this mod in other, then
            # no refinement
            if not mod_found:
                return False
        # Make sure that the relationships and activities match
        # TODO: Develop an activity hierarchy? In which kinaseactivity is a
        # refinement of activity, for example.
        if self.activity == other.activity and \
           self.relationship == other.relationship:
               return True
        else:
            return False

    def __str__(self):
        s = ("ActivityModification(%s, %s, %s, %s, %s)" %
                (self.monomer.name, self.mod, self.mod_pos, self.relationship,
                 self.activity))
        return s


class ActivatingSubstitution(Statement):
    """Statement representing the activation of a protein as a result
    of a residue substitution"""

    def __init__(self, monomer, wt_residue, pos, sub_residue, activity, rel,
                 evidence=None):
        super(ActivatingSubstitution, self).__init__(evidence)
        self.monomer = monomer
        self.wt_residue = wt_residue
        self.pos = pos
        self.sub_residue = sub_residue
        self.activity = activity
        self.rel = rel

    def matches_key(self):
        return (type(self), self.monomer.matches_key(), self.wt_residue,
                self.pos, self.sub_residue, self.activity)

    def agent_list(self):
        return [self.monomer]

    def monomers_interactions_only(self, agent_set):
        pass

    def assemble_interactions_only(self, model, agent_set):
        pass

    def __str__(self):
        s = ("ActivatingSubstitution(%s, %s, %s, %s, %s, %s)" %
                (self.monomer.name, self.wt_residue, self.pos,
                 self.sub_residue, self.activity, self.rel))
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
        return (type(self), self.gef.matches_key(), self.gef_activity,
                self.ras.matches_key())

    def agent_list(self):
        return [self.gef, self.ras]

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
        return (type(self), self.gap.matches_key(), self.gap_activity,
                self.ras.matches_key())

    def agent_list(self):
        return [self.gap, self.ras]

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
        return (type(self), tuple(m.matches_key() for m in self.members))

    def agent_list(self):
        return self.members

    def __str__(self):
        s = ("Complex(%s)" % ([m.name for m in self.members]))
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
                if self_agent.refinement_of(other_agent, eh, mh):
                    self_match_indices.add(self_agent_ix)
                    break
        if len(self_match_indices) != len(other.members):
            return False
        else:
            return True

