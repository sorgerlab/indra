from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import re
import rdflib
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
    patterns = ['http://www.openbel.org/bel/[pragm]_([A-Za-z]+)_.*',
                'http://www.openbel.org/bel/[a-z]+_[pr]_([A-Za-z]+)_.*',
                'http://www.openbel.org/bel/[a-z]+_complex_([A-Za-z]+)_.*',
                'http://www.openbel.org/bel/complex_([A-Za-z]+)_.*']
    for pr in patterns:
        match = re.match(pr, uri)
        if match is not None:
            return match.groups()[0]
    return None

def term_from_uri(uri):
    """Removes prepended URI information from terms."""
    if uri is None:
        return None
    # This insures that if we get a Literal with an integer value (as we
    # do for modification positions), it will get converted to a string,
    # not an integer.
    if isinstance(uri, rdflib.Literal):
        uri = str(uri.toPython())
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

class BelRdfProcessor(object):
    """The BelRdfProcessor extracts INDRA Statements from a BEL RDF model.

    Parameters
    ----------
    g : rdflib.Graph
        An RDF graph object containing the BEL model.

    Attributes
    ----------
    g : rdflib.Graph
        An RDF graph object containing the BEL model.
    statements : list[indra.statements.Statement]
        A list of extracted INDRA Statements representing direct mechanisms.
        This list should be used for assembly in INDRA.
    indirect_stmts : list[indra.statements.Statement]
        A list of extracted INDRA Statements representing indirect mechanisms.
        This list should be used for assembly or model checking in INDRA.
    converted_direct_stmts : list[str]
        A list of all direct BEL statements, as strings, that were converted
        into INDRA Statements.
    converted_indirect_stmts : list[str]
        A list of all indirect BEL statements, as strings, that were converted
        into INDRA Statements.
    degenerate_stmts : list[str]
        A list of degenerate BEL statements, as strings, in the BEL model.
    all_direct_stmts : list[str]
        A list of all BEL statements representing direct interactions,
        as strings, in the BEL model.
    all_indirect_stmts : list[str]
        A list of all BEL statements that represent indirect interactions,
        as strings, in the BEL model.
    """
    def __init__(self, g):
        self.g = g
        self.statements = []
        self.indirect_stmts = []
        self.converted_direct_stmts = []
        self.converted_indirect_stmts = []
        self.degenerate_stmts = []
        self.all_direct_stmts = []
        self.all_indirect_stmts = []

    def get_modifications(self):
        """Extract INDRA Modification Statements from BEL.

        Two SPARQL patterns are used for extracting Modifications from BEL:

        - q_phospho1 assumes that the subject is an AbundanceActivity, which
          increases/decreases a ModifiedProteinAbundance.

          Examples:

              kinaseActivity(proteinAbundance(HGNC:IKBKE))
              directlyIncreases
              proteinAbundance(HGNC:IRF3,proteinModification(P,S,385))

              phosphataseActivity(proteinAbundance(HGNC:DUSP4))
              directlyDecreases
              proteinAbundance(HGNC:MAPK1,proteinModification(P,T,185))

        - q_phospho2 assumes that the subject is a ProteinAbundance which
          increases/decreases a ModifiedProteinAbundance.

          Examples:

              proteinAbundance(HGNC:NGF) increases
              proteinAbundance(HGNC:NFKBIA,proteinModification(P,Y,42))

              proteinAbundance(HGNC:FGF1) decreases
              proteinAbundance(HGNC:RB1,proteinModification(P))
        """

        # Get statements where the subject is an activity
        q_phospho1 = prefixes + """
            SELECT ?enzName ?substrateName ?mod ?pos
                   ?stmt ?enzyme ?substrate ?rel
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?subject a belvoc:AbundanceActivity .
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
        # Get statements where the subject is a protein abundance
        q_phospho2 = prefixes + """
            SELECT ?enzName ?substrateName ?mod ?pos
                   ?stmt ?enzyme ?substrate ?rel
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?enzyme .
                ?stmt belvoc:hasObject ?object .
                ?enzyme a belvoc:ProteinAbundance .
                ?enzyme belvoc:hasConcept ?enzName .
                ?object a belvoc:ModifiedProteinAbundance .
                ?object belvoc:hasModificationType ?mod .
                ?object belvoc:hasChild ?substrate .
                ?substrate belvoc:hasConcept ?substrateName .
                OPTIONAL { ?object belvoc:hasModificationPosition ?pos . }
            }
        """
        for q_phospho in (q_phospho1, q_phospho2):
            # Run the query
            res_phospho = self.g.query(q_phospho)

            for stmt in res_phospho:
                # Parse out the elements of the query
                evidence = self._get_evidence(stmt[4])
                enz = self._get_agent(stmt[0], stmt[5])
                #act_type = name_from_uri(stmt[1])
                sub = self._get_agent(stmt[1], stmt[6])
                mod = term_from_uri(stmt[2])
                residue = self._get_residue(mod)
                mod_pos = term_from_uri(stmt[3])
                stmt_str = strip_statement(stmt[4])
                # Get the relationship (increases/decreases, etc.)
                rel = term_from_uri(stmt[7])
                if rel == 'DirectlyIncreases' or rel == 'DirectlyDecreases':
                    is_direct = True
                else:
                    is_direct = False

                # Build the INDRA statement
                # Handle PhosphorylationSerine, etc.
                if mod.startswith('Phosphorylation'):
                    modtype = 'phosphorylation'
                else:
                    modtype = mod.lower()
                # Get the class and invert if needed
                modclass = modtype_to_modclass[modtype]
                if rel == 'DirectlyDecreases' or rel == 'Decreases':
                    modclass = modclass_to_inverse[modclass]
                stmt = modclass(enz, sub, residue, mod_pos, evidence)
                if is_direct:
                    self.statements.append(stmt)
                    self.converted_direct_stmts.append(stmt_str)
                else:
                    self.converted_indirect_stmts.append(stmt_str)
                    self.indirect_stmts.append(stmt)
        return

    def get_composite_activating_mods(self):
        """Extract INDRA ActiveForm Statements with multiple mods from BEL.

        The SPARQL pattern used for extraction from BEL looks for a
        CompositeAbundance as subject where two constituents of the composite
        are both ModifiedProteinAbundances. The object has to be a
        Activity of a ProteinAbundance.

        Examples:

            compositeAbundance(
            proteinAbundance(PFH:"AKT Family",proteinModification(P,S,473)),
            proteinAbundance(PFH:"AKT Family",proteinModification(P,T,308)))
            directlyIncreases
            kinaseActivity(proteinAbundance(PFH:"AKT Family"))
        """
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
            evidence = self._get_evidence(stmt[7])
            # Parse out the elements of the query
            species = self._get_agent(stmt[0], stmt[8])
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
            self.converted_direct_stmts.append(stmt_str)
            st = ActiveForm(species, act_type, is_active, evidence)
            self.statements.append(st)

    def get_activating_mods(self):
        """Extract INDRA ActiveForm Statements with a single mod from BEL.

        The SPARQL pattern used for extraction from BEL looks for a
        ModifiedProteinAbundance as subject and an Activiy of a
        ProteinAbundance as object.

        Examples:

            proteinAbundance(HGNC:INSR,proteinModification(P,Y))
            directlyIncreases
            kinaseActivity(proteinAbundance(HGNC:INSR))
        """
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
            evidence = self._get_evidence(stmt[5])
            # Parse out the elements of the query
            species = self._get_agent(stmt[0], stmt[6])
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
            self.converted_direct_stmts.append(stmt_str)
            st = ActiveForm(species, act_type, is_active, evidence)
            self.statements.append(st)

    def get_complexes(self):
        """Extract INDRA Complex Statements from BEL.

        The SPARQL query used to extract Complexes looks for ComplexAbundance
        terms and their constituents. This pattern is distinct from other
        patterns in this processor in that it queries for terms, not
        full statements.

        Examples:

            complexAbundance(proteinAbundance(HGNC:PPARG),
            proteinAbundance(HGNC:RXRA))
            decreases
            biologicalProcess(MESHPP:"Insulin Resistance")
        """
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
            ev = self._get_evidence(stmt_uri)
            for e in ev:
                e.epistemics['direct'] = True
            cmplx_name = term_from_uri(stmt[0])
            cmplx_id = stmt_uri + '#' + cmplx_name
            child = self._get_agent(stmt[1], stmt[2])
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
        """Extract INDRA ActiveForm Statements based on a mutation from BEL.

        The SPARQL pattern used to extract ActiveForms due to mutations look
        for a ProteinAbundance as a subject which has a child encoding the
        amino acid substitution. The object of the statement is an
        ActivityType of the same ProteinAbundance, which is either increased
        or decreased.

        Examples:

            proteinAbundance(HGNC:NRAS,substitution(Q,61,K))
            directlyIncreases
            gtpBoundActivity(proteinAbundance(HGNC:NRAS))

            proteinAbundance(HGNC:TP53,substitution(F,134,I))
            directlyDecreases
            transcriptionalActivity(proteinAbundance(HGNC:TP53))
        """
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
            evidence = self._get_evidence(stmt[4])
            # Parse out the elements of the query
            enz = self._get_agent(stmt[0], stmt[5])
            sub_expr = term_from_uri(stmt[1])
            act_type = term_from_uri(stmt[2]).lower()
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
            self.converted_direct_stmts.append(stmt_str)
            st = ActiveForm(enz, act_type, is_active, evidence)
            self.statements.append(st)

    def get_activation(self):
        """Extract INDRA Inhibition/Activation Statements from BEL.

        The SPARQL query used to extract Activation Statements looks for
        patterns in which the subject is is an ActivityType
        (of a ProtainAbundance) or an Abundance (of a small molecule).
        The object has to be the ActivityType (typically of a
        ProteinAbundance) which is either increased or decreased.

        Examples:

            abundance(CHEBI:gefitinib) directlyDecreases
            kinaseActivity(proteinAbundance(HGNC:EGFR))

            kinaseActivity(proteinAbundance(HGNC:MAP3K5))
            directlyIncreases kinaseActivity(proteinAbundance(HGNC:MAP2K7))

        This pattern covers the extraction of Gap/Gef and GtpActivation
        Statements, which are recognized by the object activty or the
        subject activity, respectively, being `gtpbound`.

        Examples:

            catalyticActivity(proteinAbundance(HGNC:RASA1))
            directlyDecreases
            gtpBoundActivity(proteinAbundance(PFH:"RAS Family"))

            catalyticActivity(proteinAbundance(HGNC:SOS1))
            directlyIncreases
            gtpBoundActivity(proteinAbundance(HGNC:HRAS))

            gtpBoundActivity(proteinAbundance(HGNC:HRAS))
            directlyIncreases
            catalyticActivity(proteinAbundance(HGNC:TIAM1))
        """
        q_stmts = prefixes + """
            SELECT ?subjName ?subjActType ?rel ?objName ?objActType
                   ?stmt ?subj ?obj
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?subj .
                {?subj belvoc:hasActivityType ?subjActType .
                 ?subj belvoc:hasChild ?subjProt .
                 ?subjProt belvoc:hasConcept ?subjName .}
                UNION
                {?subj a belvoc:Abundance .
                 ?subj belvoc:hasConcept ?subjName .}
                ?stmt belvoc:hasObject ?obj .
                ?obj belvoc:hasActivityType ?objActType .
                ?obj belvoc:hasChild ?objProt .
                ?objProt belvoc:hasConcept ?objName .
                FILTER (?rel = belvoc:DirectlyIncreases ||
                        ?rel = belvoc:DirectlyDecreases)
            }
        """
        res_stmts = self.g.query(q_stmts)

        for stmt in res_stmts:
            evidence = self._get_evidence(stmt[5])
            subj = self._get_agent(stmt[0], stmt[6])
            subj_activity = stmt[1]
            if subj_activity:
                subj_activity = term_from_uri(stmt[1]).lower()
                subj.activity = ActivityCondition(subj_activity, True)
            rel = term_from_uri(stmt[2])
            if rel == 'DirectlyDecreases':
                is_activation = False
            else:
                is_activation = True
            obj = self._get_agent(stmt[3], stmt[7])
            obj_activity = term_from_uri(stmt[4]).lower()
            stmt_str = strip_statement(stmt[5])
            # Mark this as a converted statement
            self.converted_direct_stmts.append(stmt_str)

            # Distinguish the case when the activator is a GTPase
            # (since this may involve unique and stereotyped mechanisms)
            if subj_activity == 'gtpbound':
                if not is_activation:
                    logger.warning('GtpActivation only handles positive '
                                   'activation.')
                    continue
                self.statements.append(
                     GtpActivation(subj, obj, obj_activity, evidence))
            # If the object is a GTPase, and the subject *increases*
            # its GtpBound activity, then the subject is a GEF
            elif obj_activity == 'gtpbound' and rel == 'DirectlyIncreases':
                self.statements.append(
                        Gef(subj, obj, evidence))
            # If the object is a GTPase, and the subject *decreases*
            # its GtpBound activity, then the subject is a GAP
            elif obj_activity == 'gtpbound' and rel == 'DirectlyDecreases':
                self.statements.append(
                        Gap(subj, obj, evidence))
            # Otherwise, create a generic Activity->Activity statement
            else:
                if rel == 'DirectlyDecreases':
                    st = Inhibition(subj, obj, obj_activity, evidence)
                else:
                    st = Activation(subj, obj, obj_activity, evidence)
                self.statements.append(st)

    def get_transcription(self):
        """Extract Increase/DecreaseAmount INDRA Statements from BEL.

        Three distinct SPARQL patterns are used to extract amount
        regulations from BEL.

        - q_tscript1 searches for a subject which is a Transcription
          ActivityType of a ProteinAbundance and an object which is
          an RNAAbundance that is either increased or decreased.

          Examples:

              transcriptionalActivity(proteinAbundance(HGNC:FOXP2))
              directlyIncreases
              rnaAbundance(HGNC:SYK)

              transcriptionalActivity(proteinAbundance(HGNC:FOXP2))
              directlyDecreases
              rnaAbundance(HGNC:CALCRL)

        - q_tscript2 searches for a subject which is a ProteinAbundance
          and an object which is an RNAAbundance. Note that this pattern
          typically exists in an indirect form (i.e. increases/decreases).

          Example:

              proteinAbundance(HGNC:MTF1) directlyIncreases
              rnaAbundance(HGNC:LCN1)

        - q_tscript3 searches for a subject which is a
          ModifiedProteinAbundance, with an object which is an RNAAbundance.
          In the BEL large corpus, this pattern is found for
          subjects which are protein families or mouse/rat proteins, and
          the predicate in an indirect increase.

          Example:

              proteinAbundance(PFR:"Akt Family",proteinModification(P))
              increases
              rnaAbundance(RGD:Cald1)
        """
        q_tscript1 = prefixes + """
            SELECT ?tfName ?targetName ?stmt ?tf ?target ?rel
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?target .
                ?subject a belvoc:AbundanceActivity .
                ?subject belvoc:hasActivityType belvoc:Transcription .
                ?subject belvoc:hasChild ?tf .
                ?tf a belvoc:ProteinAbundance .
                ?tf belvoc:hasConcept ?tfName .
                ?target a belvoc:RNAAbundance .
                ?target belvoc:hasConcept ?targetName .
            }
        """
        q_tscript2 = prefixes + """
            SELECT ?tfName ?targetName ?stmt ?tf ?target ?rel
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?tf .
                ?stmt belvoc:hasObject ?target .
                ?tf a belvoc:ProteinAbundance .
                ?tf belvoc:hasConcept ?tfName .
                ?target a belvoc:RNAAbundance .
                ?target belvoc:hasConcept ?targetName .
            }
        """
        q_tscript3 = prefixes + """
            SELECT ?tfName ?targetName ?stmt ?tf ?target ?rel ?mod ?pos
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?target .
                ?subject a belvoc:ModifiedProteinAbundance .
                ?subject belvoc:hasModificationType ?mod .
                ?subject belvoc:hasChild ?tf .
                ?tf belvoc:hasConcept ?tfName .
                ?target a belvoc:RNAAbundance .
                ?target belvoc:hasConcept ?targetName .
                OPTIONAL { ?subject belvoc:hasModificationPosition ?pos . }
            }
        """
        for q_tscript in (q_tscript1, q_tscript2, q_tscript3):
            res_tscript = self.g.query(q_tscript)
            for stmt in res_tscript:
                # Get modifications on the subject, if any
                if q_tscript == q_tscript1:
                    tf = self._get_agent(stmt[0], stmt[3])
                    tf.activity = ActivityCondition('transcription', True)
                elif q_tscript == q_tscript3:
                    mod = term_from_uri(stmt[6])
                    mod_pos = term_from_uri(stmt[7])
                    mc = self._get_mod_condition(mod, mod_pos)
                    if mc is None:
                        continue
                    tf = self._get_agent(stmt[0], stmt[3])
                    tf.mods = mods=[mc]
                else:
                    tf = self._get_agent(stmt[0], stmt[3])
                # Parse out the elements of the query
                evidence = self._get_evidence(stmt[2])
                target = self._get_agent(stmt[1], stmt[4])
                stmt_str = strip_statement(stmt[2])
                # Get the relationship (increases/decreases, etc.)
                rel = term_from_uri(stmt[5])
                if rel == 'DirectlyIncreases' or rel == 'DirectlyDecreases':
                    is_direct = True
                else:
                    is_direct = False
                # Build the INDRA statement
                stmt = None
                if rel == 'DirectlyIncreases' or rel == 'Increases':
                    stmt = IncreaseAmount(tf, target, evidence)
                elif rel == 'DirectlyDecreases' or rel == 'Decreases':
                    stmt = DecreaseAmount(tf, target, evidence)
                # If we've matched a pattern, mark this as a converted statement
                if stmt is not None:
                    if is_direct:
                        self.statements.append(stmt)
                        self.converted_direct_stmts.append(stmt_str)
                    else:
                        self.indirect_stmts.append(stmt)
                        self.converted_indirect_stmts.append(stmt_str)

    def get_conversions(self):
        """Extract Conversion INDRA Statements from BEL.


        The SPARQL query used to extract Conversions searches for
        a subject (controller) which is an AbundanceActivity
        which directlyIncreases a Reaction with a given list of
        Reactants and Products.

        Examples:

            catalyticActivity(proteinAbundance(HGNC:HMOX1))
            directlyIncreases
            reaction(reactants(abundance(CHEBI:heme)),
            products(abundance(SCHEM:Biliverdine),
            abundance(CHEBI:"carbon monoxide")))
        """
        query = prefixes + """
            SELECT DISTINCT ?controller ?controllerName ?controllerActivity
                ?product ?productName ?reactant ?reactantName ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship ?rel .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?rxn .
                ?subject a belvoc:AbundanceActivity .
                ?subject belvoc:hasActivityType ?controllerActivity .
                ?subject belvoc:hasChild ?controller .
                ?controller belvoc:hasConcept ?controllerName .
                ?rxn a belvoc:Reaction .
                ?rxn belvoc:hasChild ?reactants .
                ?reactants rdfs:label ?reactLabel .
                FILTER (regex(?reactLabel, "^reactants.*"))
                ?rxn belvoc:hasChild ?products .
                ?products rdfs:label ?prodLabel .
                FILTER (regex(?prodLabel, "^products.*"))
                ?reactants belvoc:hasChild ?reactant .
                ?products belvoc:hasChild ?product .
                ?reactant belvoc:hasConcept ?reactantName .
                ?product belvoc:hasConcept ?productName .
            }
            """
        res = self.g.query(query)
        # We need to collect all pieces of the same statement so that we can
        # collect multiple reactants and products
        stmt_map = collections.defaultdict(list)
        for stmt in res:
            stmt_map[stmt[-1]].append(stmt)
        for stmts in stmt_map.values():
            # First we get the shared part of the Statement
            stmt = stmts[0]
            subj = self._get_agent(stmt[1], stmt[0])
            evidence = self._get_evidence(stmt[-1])
            stmt_str = strip_statement(stmt[-1])
            # Now we collect the participants
            obj_from_map = {}
            obj_to_map = {}
            for stmt in stmts:
                reactant_name = stmt[6]
                product_name = stmt[4]
                if reactant_name not in obj_from_map:
                    obj_from_map[reactant_name] = \
                        self._get_agent(stmt[6], stmt[5])
                if product_name not in obj_to_map:
                    obj_to_map[product_name] = \
                        self._get_agent(stmt[4], stmt[3])
            obj_from = list(obj_from_map.values())
            obj_to = list(obj_to_map.values())
            st = Conversion(subj, obj_from, obj_to, evidence=evidence)
            # If we've matched a pattern, mark this as a converted statement
            self.statements.append(st)
            self.converted_direct_stmts.append(stmt_str)

    def get_all_direct_statements(self):
        """Get all directlyIncreases/Decreases BEL statements.

        This method stores the results of the query in self.all_direct_stmts
        as a list of strings. The SPARQL query used to find direct BEL
        statements searches for all statements whose predicate is either
        DirectyIncreases or DirectlyDecreases.
        """
        logger.info("Getting all direct statements...\n")
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
        self.all_direct_stmts = [strip_statement(stmt[0]) for stmt in res_stmts]

    def get_all_indirect_statements(self):
        """Get all indirect increases/decreases BEL statements.

        This method stores the results of the query in self.all_indirect_stmts
        as a list of strings. The SPARQL query used to find indirect BEL
        statements searches for all statements whose predicate is either
        Increases or Decreases.
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
        self.all_indirect_stmts = [strip_statement(stmt[0]) for stmt in res_stmts]

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

        if not self.all_direct_stmts:
            self.get_all_direct_statements()
        if not self.degenerate_stmts:
            self.get_degenerate_statements()
        if not self.all_indirect_stmts:
            self.get_all_indirect_statements()

        logger.info('')
        logger.info("Total indirect statements: %d" %
                     len(self.all_indirect_stmts))
        logger.info("Converted indirect statements: %d" %
                     len(self.converted_indirect_stmts))
        logger.info(">> Unhandled indirect statements: %d" %
                     (len(self.all_indirect_stmts) -
                      len(self.converted_indirect_stmts)))
        logger.info('')
        logger.info("Total direct statements: %d" % len(self.all_direct_stmts))
        logger.info("Converted direct statements: %d" %
                    len(self.converted_direct_stmts))
        logger.info("Degenerate direct statements: %d" %
                    len(self.degenerate_stmts))
        logger.info(">> Unhandled direct statements: %d" %
                     (len(self.all_direct_stmts) -
                      len(self.converted_direct_stmts) -
                      len(self.degenerate_stmts)))

        logger.info('')
        logger.info("--- Unhandled direct statements ---------")
        for stmt in self.all_direct_stmts:
            if not (stmt in self.converted_direct_stmts or
                    stmt in self.degenerate_stmts):
                logger.info(stmt)
        logger.info('')
        logger.info("--- Unhandled indirect statements ---------")
        for stmt in self.all_indirect_stmts:
            if not (stmt in self.converted_indirect_stmts or
                    stmt in self.degenerate_stmts):
                logger.info(stmt)

    def print_statements(self):
        """Print all extracted INDRA Statements."""
        logger.info('--- Direct INDRA statements ----------')
        for i, stmt in enumerate(self.statements):
            logger.info("%s: %s" % (i, stmt))
        logger.info('--- Indirect INDRA statements ----------')
        for i, stmt in enumerate(self.indirect_stmts):
            logger.info("%s: %s" % (i, stmt))

    @staticmethod
    def _get_agent(concept, entity):
        name = term_from_uri(concept)
        namespace = namespace_from_uri(entity)
        db_refs = {}
        if namespace == 'HGNC':
            agent_name = name
            hgnc_id = hgnc_client.get_hgnc_id(name)
            if hgnc_id is not None:
                db_refs['HGNC'] = str(hgnc_id)
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                if up_id:
                    db_refs['UP'] = up_id
                else:
                    logger.warning('HGNC entity %s with HGNC ID %s has no '
                                   'corresponding Uniprot ID.' %
                                   (name, hgnc_id))
            else:
                logger.warning("Couldn't get HGNC ID for HGNC symbol %s" %
                               name)
        elif namespace in ('MGI', 'RGD'):
            agent_name = name
            db_refs[namespace] = name
        elif namespace in ('PFH', 'SFAM'):
            indra_name = bel_to_indra.get(name)
            db_refs[namespace] = name
            if indra_name is None:
                agent_name = name
                msg = 'Could not find mapping for BEL family: %s' % name
                logger.warning(msg)
            else:
                db_refs['FPLX'] = indra_name
                db_refs['TEXT'] = name
                agent_name = indra_name
        elif namespace in ('NCH', 'SCOMP'):
            indra_name = bel_to_indra.get(name)
            db_refs[namespace] = name
            if indra_name is None:
                agent_name = name
                msg = 'Could not find mapping for BEL complex: %s' % name
                logger.warning(msg)
            else:
                db_refs['FPLX'] = indra_name
                db_refs['TEXT'] = name
                agent_name = indra_name
        elif namespace == 'CHEBI':
            chebi_id = chebi_name_id.get(name)
            if chebi_id:
                db_refs['CHEBI'] = chebi_id
            else:
                logger.warning('CHEBI name %s not found in map.' % name)
            agent_name = name
        elif namespace == 'EGID':
            hgnc_id = hgnc_client.get_hgnc_from_entrez(name)
            db_refs['EGID'] = name
            if hgnc_id is not None:
                db_refs['HGNC'] = str(hgnc_id)
                agent_name = hgnc_client.get_hgnc_name(hgnc_id)
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                if up_id:
                    db_refs['UP'] = up_id
                else:
                    logger.warning('HGNC entity %s with HGNC ID %s has no '
                                   'corresponding Uniprot ID.' %
                                   (name, hgnc_id))
            else:
                logger.warning('Could not map EGID%s to HGNC.' % name)
                agent_name = 'E%s' % name
        else:
            logger.warning('Unhandled entity namespace: %s' % namespace)
            print('%s, %s' % (concept, entity))
            agent_name = name
        agent = Agent(agent_name, db_refs=db_refs)
        return agent

    def _get_evidence(self, statement):
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

        # Query for directness
        q_direct = prefixes + """
            SELECT ?predicate
            WHERE {
                <%s> belvoc:hasRelationship ?predicate .
            }
        """ % statement.format()
        res_direct = self.g.query(q_direct)
        epistemics = {}
        if res_direct:
            rel = term_from_uri(list(res_direct)[0][0])
            if rel in ('DirectlyDecreases', 'DirectlyIncreases'):
                epistemics['direct'] = True
            if rel in ('Decreases', 'Increases'):
                epistemics['direct'] = False

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
                                  annotations=annotations, epistemics=epistemics)
                    evs.append(ev)
                else:
                    logger.warning('Could not parse citation: %s' % citation)
        if not evs:
            evs = [Evidence(source_api='bel', source_id=statement,
                            annotations=annotations, epistemics=epistemics)]
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
        if not mod:
            return None
        if mod.startswith('Phosphorylation'):
            mc = ModCondition('phosphorylation')
        else:
            mc = ModCondition(mod.lower())
        mc.residue = BelRdfProcessor._get_residue(mod)
        mc.position = mod_pos
        return mc


def _build_famplex_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../../resources/famplex_map.tsv')
    bel_to_indra = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        namespace = row[0]
        entry = row[1]
        indra_name = row[2]
        if namespace == 'BEL':
            bel_to_indra[entry] = indra_name
    return bel_to_indra


def _build_chebi_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../../resources/bel_chebi_map.tsv')
    chebi_name_id = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        chebi_name = row[0]
        chebi_id = row[1]
        chebi_name_id[chebi_name] = chebi_id
    return chebi_name_id

bel_to_indra = _build_famplex_map()
chebi_name_id = _build_chebi_map()
