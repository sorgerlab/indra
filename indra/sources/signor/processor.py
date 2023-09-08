"""
An input processor for the SIGNOR database: a database of causal relationships
between biological entities.

See publication:

Perfetto et al., "SIGNOR: a database of causal relationships between
biological entities," Nucleic Acids Research, Volume 44, Issue D1, 4
January 2016, Pages D548-D554. https://doi.org/10.1093/nar/gkv1048
"""
import re
import logging
from copy import deepcopy
from collections import Counter
from os.path import join, dirname
import tqdm
from indra.statements import *
from indra.util import read_unicode_csv
from indra.resources import get_resource_path
from indra.ontology.standardize import standardize_name_db_refs, \
    get_standard_agent
from indra.sources.reach.processor import parse_amino_acid_string
from indra.databases import hgnc_client, uniprot_client, chebi_client
from indra.databases.identifiers import ensure_prefix

logger = logging.getLogger(__name__)


def _read_famplex_map():
    fname = get_resource_path('famplex_map.tsv')
    raw_map = read_unicode_csv(fname, '\t')

    m = {}
    for row in raw_map:
        m[(row[0], row[1])] = row[2]
    return m


famplex_map = _read_famplex_map()


_default_csv_file = join(dirname(__file__), '..', '..', '..', 'data',
                         'all_data_23_09_17.csv')


_type_db_map = {
    ('antibody', None): None,
    ('protein', 'UNIPROT'): 'UP',
    ('complex', 'SIGNOR'): 'SIGNOR',
    ('proteinfamily', 'SIGNOR'): 'SIGNOR',
    ('smallmolecule', 'PUBCHEM'): 'PUBCHEM',
    ('pathway', None): None,
    ('phenotype', 'SIGNOR'): 'SIGNOR',
    ('stimulus', 'SIGNOR'): 'SIGNOR',
    ('chemical', 'PUBCHEM'): 'PUBCHEM',
    ('fusion protein', 'SIGNOR'): 'SIGNOR',
    ('chemical', 'ChEBI'): 'CHEBI',
    ('smallmolecule', 'ChEBI'): 'CHEBI',
    ('mirna', 'miRBase'): 'MIRBASE',
    ('antibody', 'DRUGBANK'): 'DRUGBANK',
    ('ncrna', 'RNAcentral'): 'RNACENTRAL',
}


_mechanism_map = {
    'catalytic activity': None,
    'oxidoreductase activity': None,
    'transcriptional activation': None,
    'transcriptional repression': None,
    'Farnesylation': Farnesylation,
    'gtpase-activating protein': Gap,
    'deacetylation': Deacetylation,
    'demethylation': Demethylation,
    'dephosphorylation': Dephosphorylation,
    'destabilization': DecreaseAmount,
    'guanine nucleotide exchange factor': Gef,
    'acetylation': Acetylation,
    'binding': Complex,
    'cleavage': None,
    'desumoylation': Desumoylation,
    'deubiquitination': Deubiquitination,
    'glycosylation': Glycosylation,
    'hydroxylation': Hydroxylation,
    'neddylation': None,
    'chemical activation': Activation,
    'chemical inhibition': Inhibition,
    'trimethylation': Methylation,
    'ubiquitination': Ubiquitination,
    'monoubiquitination': Ubiquitination,
    'polyubiquitination': Ubiquitination,
    'post transcriptional regulation': None,
    'relocalization': None, # TODO: Translocation,
    'small molecule catalysis': None,
    's-nitrosylation': None,
    'transcriptional regulation': None,
    'translation regulation': None,
    'tyrosination': None,
    'lipidation': None,
    'oxidation': None,
    'methylation': Methylation,
    'palmitoylation': Palmitoylation,
    'phosphorylation': Phosphorylation,
    'stabilization': IncreaseAmount,
    'sumoylation': Sumoylation,
}


_effect_map = {
    'down-regulates': Inhibition, # TODO: Need generic downregulation
    'down-regulates activity': Inhibition,
    'down-regulates quantity': DecreaseAmount,
    'down-regulates quantity by destabilization': DecreaseAmount,
    'down-regulates quantity by repression': DecreaseAmount,
    'form complex': Complex,
    'unknown': None,
    'up-regulates': Activation, # TODO: Need generic upregulation
    'up-regulates activity': Activation,
    'up-regulates quantity': IncreaseAmount,
    'up-regulates quantity by expression': IncreaseAmount,
    'up-regulates quantity by stabilization': IncreaseAmount
}


class SignorProcessor(object):
    """Processor for Signor dataset, available at http://signor.uniroma2.it.

    Parameters
    ----------
    data : iterator
        Iterator over rows of a SIGNOR CSV file.
    complex_map : dict
        A dict containing SIGNOR complexes, keyed by their IDs.

    Attributes
    ----------
    statements : list[indra.statements.Statements]
        A list of INDRA Statements extracted from the SIGNOR table.
    stats : dict
        A dictionary containing statistics about the processing, useful
        for determining any unprocessed entries and debugging.
    """
    def __init__(self, data, complex_map=None):
        self._data = data
        if complex_map is None:
            self.complex_map = {}
        else:
            self.complex_map = complex_map
        self.stats = {}

        # Process into statements
        self.statements = []

        # Keys missing from FamPlex map
        self.stats['famplex_missing'] = []

        # Counter listing the frequency of different mechanisms that are
        # not handled by the processor.
        self.stats['unhandled_mech_ctr'] = Counter()

        # List of SignorRow namedtuples
        # List of rows where no mechanism statements were generated.
        self.stats['no_mech_rows'] = []

        for idx, row in enumerate(tqdm.tqdm(self._data,
                                            desc='Processing SIGNOR rows')):
            row_stmts, no_mech = self._process_row(row)
            if row_stmts is None:
                continue
            if no_mech:
                self.stats['no_mech_rows'].append(row)
            self.statements.extend(row_stmts)

        # Counter listing the frequency of different MECHANISM types in the
        # list of no-mechanism rows.
        # No-mechanism rows by mechanism type
        no_mech_ctr = Counter([row.MECHANISM
                               for row in self.stats['no_mech_rows']])
        self.stats['no_mech_ctr'] = \
            sorted([(k, v) for k, v in no_mech_ctr.items()],
                   key=lambda x: x[1], reverse=True)

        # Add a Complex statement for each Signor complex
        for complex_id in tqdm.tqdm(sorted(self.complex_map.keys()),
                                    desc='Processing SIGNOR complexes'):
            agents = self._get_complex_agents(complex_id)
            if len(agents) < 2:
                logger.info('Skipping Complex %s with less than 2 members' %
                            complex_id)
                continue
            # If we returned with None, we skip this complex
            if not agents:
                continue
            ev = Evidence(source_api='signor', source_id=complex_id,
                          text='Inferred from SIGNOR complex %s' % complex_id)
            s = Complex(agents, evidence=[ev])
            self.statements.append(s)
        self._log_stats()

    def _log_stats(self):
        """Log statistics about the processing."""
        logger.info('Famplex mapping missing for %d families/complexes' %
                    len(Counter(self.stats['famplex_missing'])))
        logger.info('No mechanism rows: %d' % len(self.stats['no_mech_rows']))
        logger.info('Unhandled mechanism types: %d' %
                    len(self.stats['unhandled_mech_ctr']))

    def _get_agent(self, ent_name, ent_type, id, database):
        # Returns a list of agents corresponding to this id
        # (If it is a signor complex, returns an Agent object with complex
        # constituents as BoundConditions
        name = ent_name
        if database == 'SIGNOR' and id in self.complex_map:
            components = self.complex_map[id]
            agents = self._get_complex_agents(id)
            # Return the first agent with the remaining agents as a bound
            # condition
            agent = agents[0]
            agent.bound_conditions = \
                [BoundCondition(a, True) for a in agents[1:]]
            return agent
        elif ent_type == 'mirna' and id.startswith('URS'):
            db_refs = {'RNACENTRAL': id}
            return get_standard_agent(name, db_refs=db_refs)
        else:
            gnd_type = _type_db_map[(ent_type, database)]
            if gnd_type == 'UP':
                db_refs = process_uniprot_entry(id)
            # Map SIGNOR protein families to FamPlex families
            elif ent_type == 'proteinfamily':
                db_refs = {database: id}  # Keep the SIGNOR family ID in db_refs
                key = (database, id)
                # Use SIGNOR name unless we have a mapping in FamPlex
                famplex_id = famplex_map.get(key)
                if famplex_id is None:
                    logger.debug('Could not find %s in FamPlex map' % str(key))
                    self.stats['famplex_missing'].append(key[1])
                else:
                    db_refs['FPLX'] = famplex_id
            # Other possible groundings are PUBCHEM, SIGNOR, etc.
            elif gnd_type is not None:
                if database not in ('PUBCHEM', 'SIGNOR', 'ChEBI', 'miRBase',
                                    'DRUGBANK', 'RNAcentral'):
                    raise ValueError('Unexpected database %s' % database)
                if database == 'PUBCHEM' and id.startswith('CID:'):
                    # We take off the CID: prefix plus fix an issue with
                    # SIGNOR's format in which it leaves extra spaces around
                    # the ID, as in 'CID: 923'
                    id = id[4:].strip()
                # In older releases PubChem substance IDs were used with
                # ChEBI as the source, these were later changed to use
                # PUBCHEM
                elif database in {'ChEBI', 'PUBCHEM'} \
                        and id.startswith('SID:'):
                    gnd_type = 'PUBCHEM.SUBSTANCE'
                    id = id[4:].strip()
                db_refs = {gnd_type: id}
            # If no grounding, include as an untyped/ungrounded node
            else:
                name = ent_name
                db_refs = {}
            return get_standard_agent(name, db_refs=db_refs)

    def _recursively_lookup_complex(self, complex_id):
        """Looks up the constitutents of a complex. If any constituent is
        itself a complex, recursively expands until all constituents are
        not complexes."""
        assert complex_id in self.complex_map

        expanded_agent_strings = []
        expand_these_next = [complex_id]
        while len(expand_these_next) > 0:
            # Pop next element
            c = expand_these_next[0]
            expand_these_next = expand_these_next[1:]

            # If a complex, add expanding it to the end of the queue
            # If an agent string, add it to the agent string list immediately
            assert c in self.complex_map
            for s in self.complex_map[c]:
                if s in self.complex_map and s != c:
                    expand_these_next.append(s)
                else:
                    expanded_agent_strings.append(s)
        return expanded_agent_strings

    def _get_complex_agents(self, complex_id):
        """Returns a list of agents corresponding to each of the constituents
        in a SIGNOR complex."""
        agents = []
        components = self._recursively_lookup_complex(complex_id)

        for c in components:
            db_refs = {}
            if c.startswith('CHEBI'):
                db_refs['CHEBI'] = c
                name = chebi_client.get_chebi_name_from_id(c)
            else:
                if not c.startswith('SIGNOR'):
                    name = uniprot_client.get_gene_name(c, web_fallback=False)
                else:
                    name = None
                if name is None:
                    db_refs['SIGNOR'] = c
                else:
                    db_refs['UP'] = c
                    hgnc_id = uniprot_client.get_hgnc_id(c)
                    if hgnc_id:
                        name = hgnc_client.get_hgnc_name(hgnc_id)
                        db_refs['HGNC'] = hgnc_id

                famplex_key = ('SIGNOR', c)
                if famplex_key in famplex_map:
                    db_refs['FPLX'] = famplex_map[famplex_key]
                    if not name:
                        # Set agent name to Famplex name if
                        # the Uniprot name is not available
                        name = db_refs['FPLX']
                elif not name:
                    # We neither have a Uniprot nor Famplex grounding
                    logger.debug('Have neither a Uniprot nor Famplex grounding '
                                 'for "%s" in complex %s' % (c, complex_id))
                    self.stats['famplex_missing'].append(c)
                    if not name:
                        # Set the agent name to the Signor name if neither the
                        # Uniprot nor Famplex names are available
                        name = db_refs['SIGNOR']
            assert name is not None
            agents.append(Agent(name, db_refs=db_refs))
        return agents


    @staticmethod
    def _get_evidence(row):
        # Get epistemics (direct/indirect)
        epistemics = {}
        epistemics['direct'] = True if row.DIRECT == 'YES' else False
        # Get annotations
        _n = lambda s: s if s else None
        # TODO: Refactor to exclude keys that are just Nones
        annotations = {
                'SEQUENCE': _n(row.SEQUENCE),
                'MODULATOR_COMPLEX': _n(row.MODULATOR_COMPLEX),
                'TARGET_COMPLEX': _n(row.TARGET_COMPLEX),
                'MODIFICATIONA': _n(row.MODIFICATIONA),
                'MODASEQ': _n(row.MODASEQ),
                'MODIFICATIONB': _n(row.MODIFICATIONB),
                'MODBSEQ': _n(row.MODBSEQ),
                'NOTES': _n(row.NOTES),
                'ANNOTATOR': _n(row.ANNOTATOR)}
        context = BioContext()
        if row.TAX_ID and row.TAX_ID != '-1':
            context.species = get_ref_context('TAXONOMY', row.TAX_ID)
        # NOTE: do we know if this is always a cell type, or can it be
        # a cell line?
        if row.CELL_DATA:
            # FIXME: we currently can't handle multiple pieces so we take
            # the first
            entry = row.CELL_DATA.split(';')[0]
            db_name, db_id = entry.split(':')
            context.cell_type = get_ref_context(db_name, db_id)
        # NOTE: is it okay to map this to organ?
        if row.TISSUE_DATA:
            # FIXME: we currently can't handle multiple pieces so we take
            # the first
            entry = row.TISSUE_DATA.split(';')[0]
            db_name, db_id = entry.split(':')
            context.organ = get_ref_context(db_name, db_id)
        # This is so that we don't add a blank BioContext as context and rather
        # just add None
        if not context:
            context = None

        # PMID is sometimes missing and sometimes other/Other, which we
        # don't represent
        if not row.PMID or row.PMID in {'other', 'Other'}:
            pmid = None
            text_refs = {}
        # These are regular PMIDs
        elif re.match(r'(\d+)', row.PMID):
            pmid = row.PMID
            text_refs = {'PMID': pmid}
        # Sometimes we get PMC IDs
        elif row.PMID.startswith('PMC'):
            pmid = None
            text_refs = {'PMCID': row.PMID}
        # Sometimes it's an NCBI Book
        elif row.PMID.startswith('NBK'):
            pmid = None
            text_refs = {'NCBIBOOK': row.PMID}
        # We log any other suspicious unhandled IDs
        else:
            logger.info('Invalid PMID: %s' % row.PMID)
            pmid = None
            text_refs = {}
        return Evidence(source_api='signor', source_id=row.SIGNOR_ID,
                        pmid=pmid, text=row.SENTENCE,
                        text_refs=text_refs, epistemics=epistemics,
                        annotations=annotations, context=context)

    def _process_row(self, row):
        agent_a = self._get_agent(row.ENTITYA, row.TYPEA, row.IDA,
                                  row.DATABASEA)
        agent_b = self._get_agent(row.ENTITYB, row.TYPEB, row.IDB,
                                  row.DATABASEB)
        if not agent_a.name or not agent_b.name:
            return None, None

        evidence = SignorProcessor._get_evidence(row)
        stmts = []
        no_mech = False

        # First, check for EFFECT/MECHANISM pairs giving rise to a single
        # mechanism
        # Transcriptional regulation + (up or down)
        if row.MECHANISM == 'transcriptional regulation' and \
           row.EFFECT in ('up-regulates', 'up-regulates quantity',
                          'up-regulates quantity by expression',
                          'down-regulates', 'down-regulates quantity',
                          'down-regulates quantity by repression'):
            stmt_type = IncreaseAmount if row.EFFECT.startswith('up') \
                                       else DecreaseAmount
            # Since this is a transcriptional regulation, apply a
            # transcriptional activity condition to the subject
            ac = ActivityCondition('transcription', True)
            agent_a.activity = ac
            # Create the statement
            stmts.append(stmt_type(agent_a, agent_b, evidence=evidence))
        # Stabilization + up
        elif row.MECHANISM == 'stabilization' and \
             row.EFFECT in ('up-regulates', 'up-regulates quantity',
                            'up-regulates quantity by stabilization'):
            stmts.append(IncreaseAmount(agent_a, agent_b, evidence=evidence))
        # Destabilization + down
        elif row.MECHANISM == 'destabilization' and \
             row.EFFECT in ('down-regulates', 'down-regulates quantity',
                            'down-regulates quantity by destabilization'):
            stmts.append(DecreaseAmount(agent_a, agent_b, evidence=evidence))
        # Chemical activation + up
        elif row.MECHANISM == 'chemical activation' and \
             row.EFFECT in ('up-regulates', 'up-regulates activity'):
            stmts.append(Activation(agent_a, agent_b, evidence=evidence))
        # Chemical inhibition + down
        elif row.MECHANISM == 'chemical inhibition' and \
             row.EFFECT in ('down-regulates', 'down-regulates activity'):
            stmts.append(Inhibition(agent_a, agent_b, evidence=evidence))
        # Binding + Form complex
        elif row.MECHANISM == 'binding' and row.EFFECT == 'form complex':
            stmts.append(Complex([agent_a, agent_b], evidence=evidence))
        # The above mechanism/effect combinations should be the only types
        # giving rise to statements of the same type with same args.
        # They also can't give rise to any active form statements; therefore
        # we have gotten all the statements we will get and can return.
        if stmts:
            return (stmts, False)

        # If we have a different effect/mechanism combination, we can now make
        # them separately without risk of redundancy.
        # Get the effect statement type:
        effect_stmt_type = _effect_map[row.EFFECT]
        # Get the mechanism statement type.
        if row.MECHANISM:
            if row.MECHANISM not in _mechanism_map:
                logger.debug('Unhandled mechanism type: %s' % row.MECHANISM)
                self.stats['unhandled_mech_ctr'][row.MECHANISM] += 1
                mech_stmt_type = None
            else:
                mech_stmt_type = _mechanism_map[row.MECHANISM]
        else:
            mech_stmt_type = None
        # (Note that either or both effect/mech stmt types may be None at this
        # point.)
        # First, create the effect statement:
        if effect_stmt_type == Complex:
            stmts.append(effect_stmt_type([agent_a, agent_b],
                                          evidence=evidence))
        elif effect_stmt_type:
            stmts.append(effect_stmt_type(agent_a, agent_b, evidence=evidence))

        # For modifications, we create the modification statement as well as
        # the appropriate active form.
        no_mech = False
        # Utility function for getting the polarity of the active form
        def af_is_activation(stmt, row):
            assert isinstance(stmt, Modification)
            # Get polarity of modification statement
            if isinstance(stmt, RemoveModification):
                stmt_polarity = -1
            else:
                stmt_polarity = 1
            # Get polarity of the effect
            if row.EFFECT.startswith('up'):
                effect_polarity = 1
            else:
                effect_polarity = -1
            return True if stmt_polarity * effect_polarity > 0 else False

        if mech_stmt_type and issubclass(mech_stmt_type, Modification):
            if not row.RESIDUE:
                # Modification
                mod_stmt = mech_stmt_type(agent_a, agent_b, None, None,
                                          evidence=evidence)
                stmts.append(mod_stmt)
                # ActiveForm
                if effect_stmt_type:
                    af_agent = deepcopy(agent_b)
                    af_agent.mods = [mod_stmt._get_mod_condition()]
                    # TODO: Currently this turns any upregulation associated
                    # with the modification into an ActiveForm (even
                    # up/down-regulations associated with amounts). This should
                    # be updated once we have a statement type relating Agent
                    # states to effects on amounts.
                    is_activation = af_is_activation(mod_stmt, row)
                    stmts.append(ActiveForm(af_agent, 'activity', is_activation,
                                            evidence=evidence))
            else:
                # Modification
                sites = _parse_residue_positions(row.RESIDUE)
                mod_stmts = [mech_stmt_type(agent_a, agent_b, site.residue,
                                            site.position,
                                            evidence=evidence)
                             for site in sites]
                stmts.extend(mod_stmts)
                # Active Form
                if effect_stmt_type:
                    mcs = [ms._get_mod_condition() for ms in mod_stmts]
                    af_agent = deepcopy(agent_b)
                    af_agent.mods = mcs
                    # TODO: See above.
                    is_activation = af_is_activation(mod_stmts[0], row)
                    stmts.append(ActiveForm(af_agent, 'activity', is_activation,
                                            evidence=evidence))
        # For Complex statements, we create an ActiveForm with a BoundCondition.
        elif mech_stmt_type == Complex:
            # Complex
            stmts.append(mech_stmt_type([agent_a, agent_b], evidence=evidence))
            # ActiveForm
            af_agent = deepcopy(agent_b)
            af_bc_agent = deepcopy(agent_a)
            af_agent.bound_conditions = [BoundCondition(af_bc_agent, True)]
            if row.EFFECT.startswith('up'):
                stmts.append(ActiveForm(af_agent, 'activity', True,
                                        evidence=evidence))
            elif row.EFFECT.startswith('down'):
                stmts.append(ActiveForm(af_agent, 'activity', False,
                                        evidence=evidence))
        # Other mechanism statement types
        elif mech_stmt_type:
            stmts.append(mech_stmt_type(agent_a, agent_b, evidence=evidence))
        # Mechanism statement type is None--marked as skipped
        else:
            no_mech = True
        return stmts, no_mech


def _parse_residue_positions(residue_field):
    # First see if this string contains two positions
    res_strs = [rs.strip() for rs in residue_field.split(';')]
    return [parse_amino_acid_string(rp) for rp in res_strs]


def get_ref_context(db_ns, db_id):
    db_id = db_id.strip()
    if db_ns in {'BTO'}:
        db_id = ensure_prefix(db_ns, db_id)
    standard_name, db_refs = standardize_name_db_refs({db_ns: db_id})
    return RefContext(standard_name, db_refs)


def process_uniprot_entry(up_id):
    """Process a UniProt entry ID into a db_refs structure."""
    # In older versions of SIGNOR, the ID was formatted as
    # P12345_PRO_12345 or P12345-1.
    # As of 4/2023, the ID is formatted as P12345-PRO_12345 or P12345-1.
    if up_id == 'P17861_P17861-2':
        up_id = 'P17861-2'
    parts = up_id.split('-')
    if len(parts) == 1:
        return {'UP': up_id}
    elif parts[1].startswith('PRO'):
        return {'UP': parts[0], 'UPPRO': parts[1]}
    else:
        return {'UP': parts[0], 'UPISO': up_id}