"""
An input processor for the SIGNOR database: a database of causal relationships
between biological entities.

See publication:

Perfetto et al., "SIGNOR: a database of causal relationships between
biological entities," Nucleic Acids Research, Volume 44, Issue D1, 4
January 2016, Pages D548-D554. https://doi.org/10.1093/nar/gkv1048
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
from io import StringIO, BytesIO
from copy import deepcopy
from os.path import join, dirname
import csv
import logging
import requests
from collections import namedtuple, Counter, defaultdict
from indra.statements import *
from indra.util import read_unicode_csv, read_unicode_csv_fileobj
from indra.databases import hgnc_client, uniprot_client


logger = logging.getLogger('signor')


_default_csv_file = join(dirname(__file__), '..', '..', 'data',
                         'all_data_23_09_17.csv')


_signor_fields = [
    'ENTITYA',
    'TYPEA',
    'IDA',
    'DATABASEA',
    'ENTITYB',
    'TYPEB',
    'IDB',
    'DATABASEB',
    'EFFECT',
    'MECHANISM',
    'RESIDUE',
    'SEQUENCE',
    'TAX_ID',
    'CELL_DATA',
    'TISSUE_DATA',
    'MODULATOR_COMPLEX',
    'TARGET_COMPLEX',
    'MODIFICATIONA',
    'MODASEQ',
    'MODIFICATIONB',
    'MODBSEQ',
    'PMID',
    'DIRECT',
    'NOTES',
    'ANNOTATOR',
    'SENTENCE',
    'SIGNOR_ID',
]


SignorRow = namedtuple('SignorRow', _signor_fields)


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
    ('smallmolecule', 'ChEBI'): 'CHEBI',
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
    'sumoylation': Sumoylation
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
    signor_csv : str
        Path to SIGNOR CSV file. If None is given (default), the
        SignorProcessor will attempt to download the full data file
        from http://signor.uniroma2.it/downloads.php.
    delimiter : str
        Field delimiter for CSV file. Defaults to semicolon ';'.

    Attributes
    ----------
    no_mech_rows: list of SignorRow namedtuples
        List of rows where no mechanism statements were generated.
    no_mech_ctr : collections.Counter
        Counter listing the frequency of different MECHANISM types in the
        list of no-mechanism rows.
    """
    def __init__(self, signor_csv=None, delimiter='\t'):
        # Get generator over the CSV file
        if signor_csv:
            data_iter = read_unicode_csv(signor_csv, delimiter=delimiter,
                                         skiprows=1)
        # If no CSV given, download directly from web
        else:
            url = 'https://signor.uniroma2.it/download_entity.php'
            res = requests.post(url, data={'organism':'human',
                                           'format':'csv',
                                           'submit':'Download'})
            if res.status_code == 200:
                # Python 2 -- csv.reader will need bytes
                if sys.version_info[0] < 3:
                    csv_io = BytesIO(res.content)
                # Python 3 -- csv.reader needs str
                else:
                    csv_io = StringIO(res.text)
                data_iter = read_unicode_csv_fileobj(csv_io,
                                                     delimiter=delimiter,
                                                     skiprows=1)
            else:
                raise Exception('Could not download Signor data.')
        # Process into a list of SignorRow namedtuples
        # Strip off any funky \xa0 whitespace characters
        self._data = [SignorRow(*[f.strip() for f in r]) for r in data_iter]
        # Process into statements
        self.statements = []
        self.no_mech_rows = []
        for row in self._data:
            row_stmts, no_mech = self._process_row(row)
            if no_mech:
                self.no_mech_rows.append(row)
            self.statements.extend(row_stmts)
        # No-mechanism rows by mechanism type
        no_mech_ctr = Counter([row.MECHANISM for row in self.no_mech_rows])
        self.no_mech_ctr = sorted([(k, v) for k, v in no_mech_ctr.items()],
                                  key=lambda x: x[1], reverse=True)

    @staticmethod
    def _get_agent(ent_name, ent_type, id, database):
        gnd_type = _type_db_map[(ent_type, database)]
        if gnd_type == 'UP':
            up_id = id
            db_refs = {'UP': up_id}
            name = uniprot_client.get_gene_name(up_id)
            hgnc_id = hgnc_client.get_hgnc_id(name)
            if hgnc_id:
                db_refs['HGNC'] = hgnc_id
        # Other possible groundings are PUBCHEM and SIGNOR
        elif gnd_type is not None:
            if database not in ('PUBCHEM', 'SIGNOR', 'ChEBI'):
                raise ValueError('Unexpected database %s' % database)
            if database == 'ChEBI':
                database = 'CHEBI'
            db_refs = {database: id}
            name = ent_name
        # If no grounding, include as an untyped/ungrounded node
        else:
            name = ent_name
            db_refs = {}
        return Agent(name, db_refs=db_refs)

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
                'TAX_ID': _n(row.TAX_ID),
                'CELL_DATA': _n(row.CELL_DATA),
                'TISSUE_DATA': _n(row.TISSUE_DATA),
                'MODULATOR_COMPLEX': _n(row.MODULATOR_COMPLEX),
                'TARGET_COMPLEX': _n(row.TARGET_COMPLEX),
                'MODIFICATIONA': _n(row.MODIFICATIONA),
                'MODASEQ': _n(row.MODASEQ),
                'MODIFICATIONB': _n(row.MODIFICATIONB),
                'MODBSEQ': _n(row.MODBSEQ),
                'NOTES': _n(row.NOTES),
                'ANNOTATOR': _n(row.ANNOTATOR)}
        return Evidence(source_api='signor', source_id=row.SIGNOR_ID,
                        pmid=row.PMID, text=row.SENTENCE,
                        epistemics=epistemics, annotations=annotations)

    @staticmethod
    def _process_row(row):
        agent_a = SignorProcessor._get_agent(row.ENTITYA, row.TYPEA, row.IDA,
                                             row.DATABASEA)
        agent_b = SignorProcessor._get_agent(row.ENTITYB, row.TYPEB, row.IDB,
                                             row.DATABASEB)
        evidence = SignorProcessor._get_evidence(row)
        stmts = []
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
        if mech_stmt_type and issubclass(mech_stmt_type, Modification):
            if not row.RESIDUE:
                # Modification
                mod_stmt = mech_stmt_type(agent_a, agent_b, None, None,
                                          evidence=evidence)
                stmts.append(mod_stmt)
                # ActiveForm
                af_agent = deepcopy(agent_b)
                af_agent.mods = [mod_stmt._get_mod_condition()]
                # TODO: Currently this turns any upregulation associated with
                # the modification into an ActiveForm (even up/down-regulations
                # associated with amounts). This should be updated once we have
                # a statement type relating Agent states to effects on amounts.
                if row.EFFECT.startswith('up'):
                    stmts.append(ActiveForm(af_agent, 'activity', True,
                                            evidence=evidence))
                elif row.EFFECT.startswith('down'):
                    stmts.append(ActiveForm(af_agent, 'activity', False,
                                            evidence=evidence))
            else:
                # Modification
                residues = _parse_residue_positions(row.RESIDUE)
                mod_stmts = [mech_stmt_type(agent_a, agent_b, res[0], res[1],
                                            evidence=evidence)
                             for res in residues]
                stmts.extend(mod_stmts)
                # Active Form
                mcs = [ms._get_mod_condition() for ms in mod_stmts]
                af_agent = deepcopy(agent_b)
                af_agent.mods = mcs
                # TODO: See above.
                if row.EFFECT.startswith('up'):
                    stmts.append(ActiveForm(af_agent, 'activity', True,
                                            evidence=evidence))
                elif row.EFFECT.startswith('down'):
                    stmts.append(ActiveForm(af_agent, 'activity', False,
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
    def _parse_respos(respos):
        # Split off the amino acid
        res = respos[0:3]
        pos = respos[3:]
        # Get the abbreviated amino acid
        res = amino_acids_reverse.get(res.lower())
        if not res:
            logger.warning("Could not get amino acid residue for "
                           "residue/position %s" % respos)
            return (None, None)
        # If there's no position, return residue only
        if not pos:
            return (res, None)
        # Make sure the position is an integer
        try:
            int(pos)
        except ValueError:
            logger.warning("Could not get valid position for residue/position "
                           "%s" % respos)
            return (None, None)
        return (res, pos)
    return [_parse_respos(rp) for rp in res_strs]

