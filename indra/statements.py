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
        # FIXME: Check db_refs!!!
        if not (self.name == other.name and \
                set(self.mods) == set(other.mods) and \
                set(self.mod_sites) == set(other.mod_sites) and \
                self.active == other.active and \
                len(self.bound_conditions) == len(other.bound_conditions)):
            return False

        # Check for corresponding bound_condition in the other Agent
        sorted_other_bcs = sorted(other.bound_conditions,
                                  key=lambda x: x.agent.name)
        # Check the state of all the Agents that the Agents are bound to
        for bc_ix, bc in enumerate(sorted(self.bound_conditions,
                                          key=lambda x: x.agent.name)):
            if not (bc.agent.matches(sorted_other_bcs[bc_ix].agent) and \
                    bc.is_bound == sorted_other_bcs[bc_ix].is_bound):
                # A mismatch!
                return False

        # Everything checks out, the two Agents match
        return True

    def matches_key(self):
        key = self.name
        key += str(set(self.mods))
        key += str(set(self.mod_sites))
        key += str(self.active)
        key += str(len(self.bound_conditions))
        key += ''.join(['%s%s' % (bc.agent.name, bc.is_bound)
                        for bc in sorted(self.bound_conditions,
                                         key=lambda x: x.agent.name)])
        return key

    def entity_matches_key(self):
        key = self.name
        return key

    def entity_matches(self, other):
        return self.entity_matches_key() == other.entity_matches_key()

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
        ev_str = 'Evidence(%s, %s, %s, %s)' % \
                 (self.source_api, self.pmid, self.annotations, self.text)
        return ev_str

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
        key = type(self).__name__
        key += self.enz.matches_key()
        key += self.sub.matches_key()
        key += str(self.mod)
        key += str(self.mod_pos)
        return key

    def entities_match_key(self):
        key = type(self).__name__
        key += self.enz.entity_matches_key()
        key += self.sub.entity_matches_key()
        return key

    def agent_list(self):
        return [self.enz, self.sub]

    def __str__(self):
        s = ("%s(%s, %s, %s, %s, %s)" %
                  (type(self).__name__, self.enz.name, self.sub.name, self.mod,
                   self.mod_pos, self.evidence))
        return s


class SelfModification(Statement):
    """Generic statement representing the self modification of a protein"""

    def __init__(self, enz, mod, mod_pos, evidence=None):
        super(SelfModification, self).__init__(evidence)
        self.enz = enz
        self.mod = mod
        self.mod_pos = mod_pos

    def __str__(self):
        s = ("%s(%s, %s, %s, %s)" %
             (type(self).__name__, self.enz.name, self.mod, self.mod_pos,
              self.evidence))
        return s

    def matches(self, other):
        if isinstance(other, SelfModification) and \
            self.enz.matches(other.enz) and \
            self.mod == other.mod and \
            self.mod_pos == other.mod_pos:
            return True
        else:
            return False

    def entities_match_key(self):
        key = type(self).__name__
        key += self.enz.entity_matches_key()
        return key


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
        key = type(self).__name__
        key += self.subj.matches_key()
        key += str(self.subj_activity)
        key += self.obj.matches_key()
        key += str(self.obj_activity)
        return key

    def entities_match_key(self):
        key = type(self).__name__
        key += self.subj.entity_matches_key()
        key += self.obj.entity_matches_key()
        return key

    def __str__(self):
        s = ("%s(%s, %s, %s, %s, %s, %s)" %
             (type(self).__name__, self.subj.name, self.subj_activity,
              self.relationship, self.obj.name, self.obj_activity,
              self.evidence))
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
        self.mod = mod
        self.mod_pos = mod_pos
        self.relationship = relationship
        self.activity = activity

    def matches_key(self):
        key = type(self).__name__
        key += self.monomer.matches_key()
        key += str(self.mod)
        key += str(self.mod_pos)
        key += str(self.relationship)
        key += str(self.activity)
        return key

    def entities_match_key(self):
        key = type(self).__name__
        key += self.monomer.entity_matches_key()
        return key

    def __str__(self):
        s = ("ActivityModification(%s, %s, %s, %s, %s, %s)" %
                (self.monomer.name, self.mod, self.mod_pos, self.relationship,
                 self.activity, self.evidence))
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
        key = type(self).__name__
        key += self.monomer.matches_key()
        key += str(self.wt_residue)
        key += str(self.pos)
        key += str(self.sub_residue)
        key += str(self.activity)
        return key

    def entities_match_key(self):
        key = type(self).__name__
        key += self.monomer.entity_matches_key()
        return key

    def monomers_interactions_only(self, agent_set):
        pass

    def assemble_interactions_only(self, model, agent_set):
        pass

    def __str__(self):
        s = ("ActivatingSubstitution(%s, %s, %s, %s, %s, %s, %s)" %
                (self.monomer.name, self.wt_residue, self.pos,
                 self.sub_residue, self.activity, self.rel, self.evidence))
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
        key = type(self).__name__
        key += self.gef.matches_key()
        key += self.gef_activity
        key += self.ras.matches_key()
        return key

    def entities_match_key(self):
        key = type(self).__name__
        key += self.gef.entity_matches_key()
        key += self.ras.entity_matches_key()
        return key

    def __str__(self):
        s = ("RasGef(%s, %s, %s, %s)" %
                (self.gef.name, self.gef_activity, self.ras.name, self.evidence))
        return s


class RasGap(Statement):
    """Statement representing the inactivation of a GTP-bound protein
    upon Gap activity."""

    def __init__(self, gap, gap_activity, ras, evidence=None):
        super(RasGap, self).__init__(evidence)
        self.gap = gap
        self.gap_activity = gap_activity
        self.ras = ras

    def matches_key(self):
        key = type(self).__name__
        key += self.gap.matches_key()
        key += self.gap_activity
        key += self.ras.matches_key()
        return key

    def entities_match_key(self):
        key = type(self).__name__
        key += self.gap.entity_matches_key()
        key += self.ras.entity_matches_key()
        return key

    def __str__(self):
        s = ("RasGap(%s, %s, %s, %s)" %
                (self.gap.name, self.gap_activity, self.ras.name, self.evidence))
        return s


class Complex(Statement):
    """Statement representing complex formation between a set of members"""

    def __init__(self, members, evidence=None):
        super(Complex, self).__init__(evidence)
        self.members = members

    def matches_key(self):
        key = type(self).__name__
        key += ''.join([m.matches_key() for m in self.members])
        return key

    def entities_match_key(self):
        key = type(self).__name__
        key += ''.join([m.entity_matches_key() for m in self.members])
        return key

    def __str__(self):
        s = ("Complex(%s, %s)" % ([m.name for m in self.members],
                                  self.evidence))
        return s


