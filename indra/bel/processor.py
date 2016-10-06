from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import logging
import collections
from requests.utils import unquote
from indra.statements import *
from indra.databases import hgnc_client
from indra.util import read_unicode_csv

logger = logging.getLogger('bel')

prefixes = """
    PREFIX belvoc: <http://www.openbel.org/vocabulary/>
    PREFIX belsc: <http://www.openbel.org/bel/>
    PREFIX belns: <http://www.openbel.org/bel/namespace/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>"""

def namespace_from_uri(uri):
    """Return the entity namespace from the URI. Examples:
    http://www.openbel.org/bel/p_HGNC_RAF1 -> HGNC
    http://www.openbel.org/bel/p_RGD_Raf1 -> RGD
    http://www.openbel.org/bel/p_PFH_MEK1/2_Family -> PFH
    """
    patterns = ['http://www.openbel.org/bel/p_([A-Za-z]+)_.*',
                'http://www.openbel.org/bel/[a-z]+_p_([A-Za-z]+)_.*',
                'http://www.openbel.org/bel/[a-z]+_complex_([A-Za-z]+)_.*',
                'http://www.openbel.org/bel/complex_([A-Za-z]+)_.*',
                'http://www.openbel.org/bel/a_([A-Za-z]+)_.*',
                'http://www.openbel.org/bel/g_([A-Za-z]+)_.*']
    for pr in patterns:
        match = re.match(pr, uri)
        if match is not None:
            return match.groups()[0]
    return None

def term_from_uri(uri):
    """Removes prepended URI information from terms."""
    if uri is None:
        return None
    # This is to handle URIs like
    # http://www.openbel.org/bel/namespace//MAPK%20Erk1/3%20Family
    # or
    # http://www.openbel.org/bel/namespace/MAPK%20Erk1/3%20Family
    # In the current implementation, the order of the patterns
    # matters.
    patterns = ['http://www.openbel.org/bel/namespace//(.*)',
                'http://www.openbel.org/vocabulary//(.*)',
                'http://www.openbel.org/bel//(.*)',
                'http://www.openbel.org/bel/namespace/(.*)',
                'http://www.openbel.org/vocabulary/(.*)',
                'http://www.openbel.org/bel/(.*)']
    for pr in patterns:
        match = re.match(pr, uri)
        if match is not None:
            term = match.groups()[0]
            term = unquote(term)
            return term
    # If none of the patterns match then the URI is actually a simple term
    # for instance a site: "341" or a substitution: "sub(V,600,E)"
    return uri

def strip_statement(uri):
    uri = uri.replace(r'http://www.openbel.org/bel/', '')
    uri = uri.replace(r'http://www.openbel.org/vocabulary/', '')
    return uri

class BelProcessor(object):
    """The BelProcessor extracts INDRA Statements from a BEL RDF model.

    Parameters
    ----------
    g : rdflib.Graph
        An RDF graph object containing the BEL model.

    Attributes
    ----------
    g : rdflib.Graph
        An RDF graph object containing the BEL model.
    statements : list[indra.statement.Statements]
        A list of extracted INDRA Statements. This list should be used for
        assembly in INDRA.
    all_stmts : list[str]
        A list of all BEL statements, as strings, in the BEL model.
    converted_stmts : list[str]
        A list of all BEL statements, as strings, that were converted into
        INDRA Statements.
    degenerate_stmts : list[str]
        A list of degenerate BEL statements, as strings, in the BEL model.
    indirect_stmts : list[str]
        A list of all BEL statements that represent indirect interactions, 
        as strings, in the BEL model.
    """
    def __init__(self, g):
        self.g = g
        self.statements = []
        self.all_stmts = []
        self.converted_stmts = []
        self.degenerate_stmts = []
        self.indirect_stmts = []

    def get_modifications(self):
        """Extract INDRA Modification Statements from BEL."""
        q_phospho = prefixes + """
            SELECT ?enzName ?actType ?substrateName ?mod ?pos
                   ?stmt ?enzyme ?substrate
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?subject a belvoc:AbundanceActivity .
                ?subject belvoc:hasActivityType ?actType .
                ?subject belvoc:hasChild ?enzyme .
                ?enzyme a belvoc:ProteinAbundance .
                ?enzyme belvoc:hasConcept ?enzName .
                ?object a belvoc:ModifiedProteinAbundance .
                ?object belvoc:hasModificationType ?mod .
                ?object belvoc:hasChild ?substrate .
                ?substrate belvoc:hasConcept ?substrateName .
                OPTIONAL { ?object belvoc:hasModificationPosition ?pos . }
            }
        """

        # Now make the PySB for the phosphorylation
        res_phospho = self.g.query(q_phospho)

        for stmt in res_phospho:
            evidence = self.get_evidence(stmt[5])
            # Parse out the elements of the query
            enz = self.get_agent(stmt[0], stmt[6])
            act_type = term_from_uri(stmt[1])
            sub = self.get_agent(stmt[2], stmt[7])
            mod = term_from_uri(stmt[3])
            residue = self._get_residue(mod)
            mod_pos = term_from_uri(stmt[4])
            stmt_str = strip_statement(stmt[5])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)

            if act_type == 'Kinase' and mod.startswith('Phosphorylation'):
                self.statements.append(
                        Phosphorylation(enz, sub, residue, mod_pos,
                                        evidence))
            elif act_type == 'Catalytic':
                if mod == 'Hydroxylation':
                    self.statements.append(
                            Hydroxylation(enz, sub, residue, mod_pos,
                                          evidence))
                elif mod == 'Sumoylation':
                    self.statements.append(
                            Sumoylation(enz, sub, residue, mod_pos,
                                        evidence))
                elif mod == 'Acetylation':
                    self.statements.append(
                            Acetylation(enz, sub, residue, mod_pos,
                                        evidence))
                elif mod == 'Ubiquitination':
                    self.statements.append(
                            Ubiquitination(enz, sub, residue, mod_pos,
                                           evidence))
                else:
                    logger.warning("Unknown modification type!")
                    logger.warning("Activity: %s, Mod: %s, Mod_Pos: %s" %
                                   (act_type, mod, mod_pos))
            else:
                logger.warning("Unknown modification type!")
                logger.warning("Activity: %s, Mod: %s, Mod_Pos: %s" %
                               (act_type, mod, mod_pos))

    def get_dephosphorylations(self):
        """Extract INDRA Dephosphorylation Statements from BEL."""
        q_phospho = prefixes + """
            SELECT ?phosName ?substrateName ?mod ?pos ?stmt
                   ?phosphatase ?substrate
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship belvoc:DirectlyDecreases .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?subject belvoc:hasActivityType belvoc:Phosphatase .
                ?subject belvoc:hasChild ?phosphatase .
                ?phosphatase a belvoc:ProteinAbundance .
                ?phosphatase belvoc:hasConcept ?phosName .
                ?object a belvoc:ModifiedProteinAbundance .
                ?object belvoc:hasModificationType ?mod .
                ?object belvoc:hasChild ?substrate .
                ?substrate belvoc:hasConcept ?substrateName .
                OPTIONAL { ?object belvoc:hasModificationPosition ?pos . }
            }
        """

        # Now make the PySB for the phosphorylation
        res_phospho = self.g.query(q_phospho)

        for stmt in res_phospho:
            evidence = self.get_evidence(stmt[4])
            # Parse out the elements of the query
            phos = self.get_agent(stmt[0], stmt[5])
            sub = self.get_agent(stmt[1], stmt[6])
            mod = term_from_uri(stmt[2])
            residue = self._get_residue(mod)
            mod_pos = term_from_uri(stmt[3])
            stmt_str = strip_statement(stmt[4])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.statements.append(
                    Dephosphorylation(phos, sub, residue, mod_pos, evidence))

    def get_composite_activating_mods(self):
        """Extract INDRA ActiveForm Statements with multiple mods from BEL."""
        # To eliminate multiple matches, we use pos1 < pos2 but this will
        # only work if the pos is given, otherwise multiple matches of
        # the same mod combination may appear in the result
        q_mods = prefixes + """
            SELECT ?speciesName ?actType ?mod1 ?pos1 ?mod2 ?pos2 ?rel ?stmt
                   ?species
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?object belvoc:hasActivityType ?actType .
                ?object belvoc:hasChild ?species .
                ?species a belvoc:ProteinAbundance .
                ?species belvoc:hasConcept ?speciesName .
                ?subject a belvoc:CompositeAbundance .
                ?subject belvoc:hasChild ?subject1 .
                ?subject1 a belvoc:ModifiedProteinAbundance .
                ?subject1 belvoc:hasModificationType ?mod1 .
                ?subject1 belvoc:hasChild ?species .
                ?subject belvoc:hasChild ?subject2 .
                ?subject2 a belvoc:ModifiedProteinAbundance .
                ?subject2 belvoc:hasModificationType ?mod2 .
                ?subject2 belvoc:hasChild ?species .
                OPTIONAL { ?subject1 belvoc:hasModificationPosition ?pos1 . }
                OPTIONAL { ?subject2 belvoc:hasModificationPosition ?pos2 . }
                FILTER ((?rel = belvoc:DirectlyIncreases ||
                        ?rel = belvoc:DirectlyDecreases) &&
                        ?pos1 < ?pos2)
            }
        """

        # Now make the PySB for the phosphorylation
        res_mods = self.g.query(q_mods)

        for stmt in res_mods:
            evidence = self.get_evidence(stmt[7])
            # Parse out the elements of the query
            species = self.get_agent(stmt[0], stmt[8])
            act_type = term_from_uri(stmt[1]).lower()
            mod1 = term_from_uri(stmt[2])
            mod_pos1 = term_from_uri(stmt[3])
            mc1 = self._get_mod_condition(mod1, mod_pos1)
            mod2 = term_from_uri(stmt[4])
            mod_pos2 = term_from_uri(stmt[5])
            mc2 = self._get_mod_condition(mod2, mod_pos2)
            species.mods = [mc1, mc2]
            rel = term_from_uri(stmt[6])
            if rel == 'DirectlyDecreases':
                is_active = False
            else:
                is_active = True
            stmt_str = strip_statement(stmt[7])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.statements.append(
                    ActiveForm(species, act_type, is_active, evidence))

    def get_activating_mods(self):
        """Extract INDRA ActiveForm Statements with a single mod from BEL."""
        q_mods = prefixes + """
            SELECT ?speciesName ?actType ?mod ?pos ?rel ?stmt ?species
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?object belvoc:hasActivityType ?actType .
                ?object belvoc:hasChild ?species .
                ?species a belvoc:ProteinAbundance .
                ?species belvoc:hasConcept ?speciesName .
                ?subject a belvoc:ModifiedProteinAbundance .
                ?subject belvoc:hasModificationType ?mod .
                ?subject belvoc:hasChild ?species .
                OPTIONAL { ?subject belvoc:hasModificationPosition ?pos . }
                FILTER (?rel = belvoc:DirectlyIncreases ||
                        ?rel = belvoc:DirectlyDecreases)
            }
        """

        # Now make the PySB for the phosphorylation
        res_mods = self.g.query(q_mods)

        for stmt in res_mods:
            evidence = self.get_evidence(stmt[5])
            # Parse out the elements of the query
            species = self.get_agent(stmt[0], stmt[6])
            act_type = term_from_uri(stmt[1]).lower()
            mod = term_from_uri(stmt[2])
            mod_pos = term_from_uri(stmt[3])
            mc = self._get_mod_condition(mod, mod_pos)
            species.mods = [mc]
            rel = term_from_uri(stmt[4])
            if rel == 'DirectlyDecreases':
                is_active = False
            else:
                is_active = True
            stmt_str = strip_statement(stmt[5])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.statements.append(
                    ActiveForm(species, act_type, is_active, evidence))

    def get_complexes(self):
        """Extract INDRA Complex Statements from BEL."""
        q_cmplx = prefixes + """
            SELECT ?complexTerm ?childName ?child ?stmt
            WHERE {
                {
                {?stmt belvoc:hasSubject ?complexTerm}
                UNION
                {?stmt belvoc:hasObject ?complexTerm .}
                UNION
                {?stmt belvoc:hasSubject ?term .
                ?term belvoc:hasChild ?complexTerm .}
                UNION
                {?stmt belvoc:hasObject ?term .
                ?term belvoc:hasChild ?complexTerm .}
                }
                ?complexTerm a belvoc:Term .
                ?complexTerm a belvoc:ComplexAbundance .
                ?complexTerm belvoc:hasChild ?child .
                ?child belvoc:hasConcept ?childName .
            }
        """
        # Run the query
        res_cmplx = self.g.query(q_cmplx)

        # Store the members of each complex in a dict of lists, keyed by the
        # term for the complex
        cmplx_dict = collections.defaultdict(list)
        cmplx_ev = {}
        for stmt in res_cmplx:
            stmt_uri = stmt[3]
            ev = self.get_evidence(stmt_uri)
            cmplx_name = term_from_uri(stmt[0])
            cmplx_id = stmt_uri + '#' + cmplx_name
            child = self.get_agent(stmt[1], stmt[2])
            cmplx_dict[cmplx_id].append(child)
            # This might be written multiple times but with the same
            # evidence
            cmplx_ev[cmplx_id] = ev
        # Now iterate over the stored complex information and create binding
        # statements
        for cmplx_id, cmplx_list in cmplx_dict.items():
            if len(cmplx_list) < 2:
                msg = 'Complex %s has less than 2 members! Skipping.' % \
                       cmplx_name
                logger.warning(msg)
            else:
                self.statements.append(Complex(cmplx_list,
                                               evidence=cmplx_ev[cmplx_id]))

    def get_activating_subs(self):
        """Extract INDRA ActiveForm Statements based on a mutation from BEL."""
        #p_HGNC_NRAS_sub_Q_61_K_DirectlyIncreases_gtp_p_HGNC_NRAS
        #p_HGNC_KRAS_sub_G_12_R_DirectlyIncreases_gtp_p_PFH_RAS_Family
        #p_HGNC_BRAF_sub_V_600_E_DirectlyIncreases_kin_p_HGNC_BRAF
        q_mods = prefixes + """
            SELECT ?enzyme_name ?sub_label ?act_type ?rel ?stmt ?subject
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?subject a belvoc:ProteinAbundance .
                ?subject belvoc:hasConcept ?enzyme_name .
                ?subject belvoc:hasChild ?sub_expr .
                ?sub_expr rdfs:label ?sub_label .
                ?object a belvoc:AbundanceActivity .
                ?object belvoc:hasActivityType ?act_type .
                ?object belvoc:hasChild ?enzyme .
                ?enzyme a belvoc:ProteinAbundance .
                ?enzyme belvoc:hasConcept ?enzyme_name .
            }
        """

        # Now make the PySB for the phosphorylation
        res_mods = self.g.query(q_mods)

        for stmt in res_mods:
            evidence = self.get_evidence(stmt[4])
            # Parse out the elements of the query
            enz = self.get_agent(stmt[0], stmt[5])
            sub_expr = term_from_uri(stmt[1])
            act_type = term_from_uri(stmt[2])
            # Parse the WT and substituted residues from the node label.
            # Strangely, the RDF for substituted residue doesn't break the
            # terms of the BEL expression down into their meaning, as happens
            # for modified protein abundances. Instead, the substitution
            # just comes back as a string, e.g., "sub(V,600,E)". This code
            # parses the arguments back out using a regular expression.
            match = re.match('sub\(([A-Z]),([0-9]*),([A-Z])\)', sub_expr)
            if match:
                matches = match.groups()
                wt_residue = matches[0]
                position = matches[1]
                sub_residue = matches[2]
            else:
                logger.warning("Could not parse substitution expression %s" %
                               sub_expr)
                continue
            mc = MutCondition(position, wt_residue, sub_residue)
            enz.mutations = [mc]
            rel = strip_statement(stmt[3])
            if rel == 'DirectlyDecreases':
                is_active = False
            else:
                is_active = True

            stmt_str = strip_statement(stmt[4])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.statements.append(
                    ActiveForm(enz, act_type, is_active, evidence))

    def get_activation(self):
        """Extract INDRA Activation Statements from BEL."""
        # Query for all statements where the activity of one protein
        # directlyIncreases the activity of another protein, without reference
        # to a modification.
        q_stmts = prefixes + """
            SELECT ?subjName ?subjActType ?rel ?objName ?objActType
                   ?stmt ?subj ?obj
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasSubject ?subj .
                ?stmt belvoc:hasObject ?obj .
                ?stmt belvoc:hasRelationship ?rel .
                ?subj belvoc:hasActivityType ?subjActType .
                ?subj belvoc:hasChild ?subjProt .
                ?subjProt belvoc:hasConcept ?subjName .
                ?obj belvoc:hasActivityType ?objActType .
                ?obj belvoc:hasChild ?objProt .
                ?objProt belvoc:hasConcept ?objName .
                FILTER (?rel = belvoc:DirectlyIncreases ||
                        ?rel = belvoc:DirectlyDecreases)
            }
        """
        res_stmts = self.g.query(q_stmts)

        for stmt in res_stmts:
            evidence = self.get_evidence(stmt[5])
            subj = self.get_agent(stmt[0], stmt[6])
            subj_activity = term_from_uri(stmt[1]).lower()
            rel = term_from_uri(stmt[2])
            if rel == 'DirectlyDecreases':
                is_activation = False
            else:
                is_activation = True
            obj = self.get_agent(stmt[3], stmt[7])
            obj_activity = term_from_uri(stmt[4]).lower()
            stmt_str = strip_statement(stmt[5])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)

            # Distinguish the case when the activator is a RasGTPase
            # (since this may involve unique and stereotyped mechanisms)
            if subj_activity == 'gtpbound':
                self.statements.append(
                     RasGtpActivation(subj, subj_activity,
                                      obj, obj_activity, is_activation,
                                      evidence))
            # If the object is a Ras-like GTPase, and the subject *increases*
            # its GtpBound activity, then the subject is a RasGEF
            elif obj_activity == 'gtpbound' and \
                 rel == 'increases':
                self.statements.append(
                        RasGef(subj, subj_activity, obj, evidence))
            # If the object is a Ras-like GTPase, and the subject *decreases*
            # its GtpBound activity, then the subject is a RasGAP
            elif obj_activity == 'gtpbound' and \
                 rel == 'decreases':
                self.statements.append(
                        RasGap(subj, subj_activity, obj, evidence))
            # Otherwise, create a generic Activity->Activity statement
            else:
                self.statements.append(
                     Activation(subj, subj_activity,
                                obj, obj_activity, is_activation, evidence))

            """
            #print "--------------------------------"
            print stmt_str
            print("This statement says that:")
            print("%s activity increases activity of %s" %
                  (subj_name, obj_name))
            print "It doesn't specify the site."
            act_mods = []
            for bps in self.statements:
                if type(bps) == ActivatingModification and \
                   bps.monomer_name == obj_name:
                    act_mods.append(bps)
            # If we know about an activation modification...
            if act_mods:
                print "However, I happen to know about the following"
                print "activating modifications for %s:" % obj_name
                for act_mod in act_mods:
                    print "    %s at %s" % (act_mod.mod, act_mod.mod_pos)
        """

    def get_all_direct_statements(self):
        """Get all directlyIncreases/Decreases BEL statements.

        Stores the results of the query in self.all_stmts.
        """
        logger.info("Getting all direct statements...\n")
        q_stmts = prefixes + """
            SELECT ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasSubject ?subj .
                ?stmt belvoc:hasObject ?obj .
                {
                  { ?subj a belvoc:AbundanceActivity . }
                  UNION
                  { ?subj a belvoc:ComplexAbundance . }
                  UNION
                  { ?subj a belvoc:ProteinAbundance . }
                  UNION
                  { ?subj a belvoc:ModifiedProteinAbundance . }
                }
                {
                  { ?obj a belvoc:AbundanceActivity . }
                  UNION
                  { ?obj a belvoc:ComplexAbundance . }
                  UNION
                  { ?obj a belvoc:ProteinAbundance . }
                  UNION
                  { ?obj a belvoc:ModifiedProteinAbundance . }
                }

                {
                  { ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases . }
                  UNION
                  { ?stmt belvoc:hasRelationship belvoc:DirectlyDecreases . }
                }
            }
        """
        q_stmts = prefixes + """
            SELECT ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                {
                  { ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases . }
                  UNION
                  { ?stmt belvoc:hasRelationship belvoc:DirectlyDecreases . }
                }
            }
        """

        res_stmts = self.g.query(q_stmts)
        self.all_stmts = [strip_statement(stmt[0]) for stmt in res_stmts]

    def get_indirect_statements(self):
        """Get all indirect increases/decreases BEL statements.

        Stores the results of the query in self.indirect_stmts.
        """
        q_stmts = prefixes + """
            SELECT ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                {
                  { ?stmt belvoc:hasRelationship belvoc:Increases . }
                  UNION
                  { ?stmt belvoc:hasRelationship belvoc:Decreases . }
                }
            }
        """

        res_stmts = self.g.query(q_stmts)
        self.indirect_stmts = [strip_statement(stmt[0]) for stmt in res_stmts]

    def get_degenerate_statements(self):
        """Get all degenerate BEL statements.

        Stores the results of the query in self.degenerate_stmts.
        """
        logger.info("Checking for 'degenerate' statements...\n")
        # Get rules of type protein X -> activity Y
        q_stmts = prefixes + """
            SELECT ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasSubject ?subj .
                ?stmt belvoc:hasObject ?obj .
                {
                  { ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases . }
                  UNION
                  { ?stmt belvoc:hasRelationship belvoc:DirectlyDecreases . }
                }
                {
                  { ?subj a belvoc:ProteinAbundance . }
                  UNION
                  { ?subj a belvoc:ModifiedProteinAbundance . }
                }
                ?subj belvoc:hasConcept ?xName .
                {
                  {
                    ?obj a belvoc:ProteinAbundance .
                    ?obj belvoc:hasConcept ?yName .
                  }
                  UNION
                  {
                    ?obj a belvoc:ModifiedProteinAbundance .
                    ?obj belvoc:hasChild ?proteinY .
                    ?proteinY belvoc:hasConcept ?yName .
                  }
                  UNION
                  {
                    ?obj a belvoc:AbundanceActivity .
                    ?obj belvoc:hasChild ?objChild .
                    ?objChild a belvoc:ProteinAbundance .
                    ?objChild belvoc:hasConcept ?yName .
                  }
                }
                FILTER (?xName != ?yName)
            }
        """
        res_stmts = self.g.query(q_stmts)

        logger.info("Protein -> Protein/Activity statements:")
        logger.info("---------------------------------------")
        for stmt in res_stmts:
            stmt_str = strip_statement(stmt[0])
            logger.info(stmt_str)
            self.degenerate_stmts.append(stmt_str)

    def print_statement_coverage(self):
        """Display how many of the direct statements have been converted.

        Also prints how many are considered 'degenerate' and not converted."""

        if not self.all_stmts:
            self.get_all_direct_statements()
        if not self.degenerate_stmts:
            self.get_degenerate_statements()
        if not self.indirect_stmts:
            self.get_indirect_statements()

        logger.info('')
        logger.info("Total indirect statements: %d" %
                     len(self.indirect_stmts))
        logger.info("Total direct statements: %d" % len(self.all_stmts))
        logger.info("Converted statements: %d" % len(self.converted_stmts))
        logger.info("Degenerate statements: %d" % len(self.degenerate_stmts))
        logger.info(">> Total unhandled statements: %d" %
                     (len(self.all_stmts) - len(self.converted_stmts) -
                     len(self.degenerate_stmts)))

        logger.info('')
        logger.info("--- Unhandled statements ---------")
        for stmt in self.all_stmts:
            if not (stmt in self.converted_stmts or
                    stmt in self.degenerate_stmts):
                logger.info(stmt)

    def print_statements(self):
        """Print all extracted INDRA Statements."""
        for i, stmt in enumerate(self.statements):
            logger.info("%s: %s" % (i, stmt))

    @staticmethod
    def get_agent(concept, entity):
        name = term_from_uri(concept)
        namespace = namespace_from_uri(entity)
        db_refs = {}
        if namespace == 'HGNC':
            agent_name = name
            hgnc_id = hgnc_client.get_hgnc_id(name)
            if hgnc_id is not None:
                db_refs['HGNC'] = str(hgnc_id)
        elif namespace in ('MGI', 'RGD'):
            agent_name = name
        elif namespace in ('PFH', 'SFAM'):
            indra_name = bel_to_indra.get(name)
            if indra_name is None:
                agent_name = name
                msg = 'Could not find mapping for BEL family: %s' % name
                logger.warning(msg)
            else:
                db_refs['BE'] = indra_name
                db_refs['TEXT'] = name
                agent_name = indra_name
        elif namespace == 'CHEBI':
            chebi_id = chebi_name_id.get(name)
            if chebi_id:
                db_refs['CHEBI'] = chebi_id
            else:
                logger.warning('CHEBI name %s not found in map.' % chebi_id)
            agent_name = name
        elif namespace == 'EGID':
            hgnc_id = hgnc_client.get_hgnc_from_entrez(name)
            if hgnc_id is not None:
                db_refs['HGNC'] = str(hgnc_id)
                agent_name = hgnc_client.get_hgnc_name(hgnc_id)
            else:
                logger.warning('Could not map EGID%s to HGNC.' % name)
                agent_name = 'E%s' % name
        else:
            logger.warning('Unhandled entity namespace: %s' % namespace)
            print('%s, %s' % (concept, entity))
            agent_name = name
        agent = Agent(agent_name, db_refs=db_refs)
        return agent

    def get_evidence(self, statement):
        evidence = None
        citation = None
        annotations = []

        # Query for all annotations of the statement
        q_annotations = prefixes + """
            SELECT ?annotation
            WHERE {
                <%s> belvoc:hasEvidence ?evidence .
                ?evidence belvoc:hasAnnotation ?annotation .
            }
        """ % statement.format()
        res_annotations = self.g.query(q_annotations)
        for stmt in res_annotations:
            annotations.append(stmt[0].format())

        # Query for evidence text and citation
        q_evidence = prefixes + """
            SELECT ?evidenceText ?citation
            WHERE {
                <%s> belvoc:hasEvidence ?evidence .
                ?evidence belvoc:hasEvidenceText ?evidenceText .
                ?evidence belvoc:hasCitation ?citation .
            }
        """ % statement.format()
        res_evidence = self.g.query(q_evidence)
        evs = []
        for stmt in res_evidence:
            text = stmt[0].toPython()
            citation = stmt[1].toPython()
            if citation is not None:
                m = re.match('.*pubmed:([0-9]+)', citation)
                if m is not None:
                    citation = m.groups()[0]
                    ev = Evidence(source_api='bel', source_id=statement,
                                  pmid=citation, text=text,
                                  annotations=annotations)
                    evs.append(ev)
                else:
                    logger.warning('Could not parse citation: %s' % citation)
        if not evs:
            evs = [Evidence(source_api='bel', source_id=statement,
                            annotations=annotations)]
        return evs

    @staticmethod
    def _get_residue(mod):
        if mod.startswith('Phosphorylation'):
            if mod == 'Phosphorylation':
                residue = None
            else:
                residue = mod[15:].lower()
                residue = get_valid_residue(residue)
        else:
            residue = None
        return residue

    @staticmethod
    def _get_mod_condition(mod, mod_pos):
        if mod.startswith('Phosphorylation'):
            mc = ModCondition('phosphorylation')
        else:
            mc = ModCondition(mod)
        mc.residue = BelProcessor._get_residue(mod)
        mc.position = mod_pos
        return mc

def _build_bel_indra_map():
    fname = os.path.dirname(os.path.abspath(__file__)) +\
                '/../resources/bel_indra_map.tsv'
    bel_to_indra = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        bel_name = row[0]
        indra_name = row[1]
        bel_to_indra[bel_name] = indra_name
    return bel_to_indra

def _build_chebi_map():
    fname = os.path.dirname(os.path.abspath(__file__)) +\
                '/../resources/bel_chebi_map.tsv'
    chebi_name_id = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        chebi_name = row[0]
        chebi_id = row[1]
        chebi_name_id[chebi_name] = chebi_id
    return chebi_name_id

bel_to_indra = _build_bel_indra_map()
chebi_name_id = _build_chebi_map()
