from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from io import StringIO
from copy import deepcopy
from os.path import join, dirname
import csv
import logging
import requests
from collections import namedtuple, Counter, defaultdict

from indra.statements import *
from indra.util import read_unicode_csv
from indra.databases import hgnc_client, uniprot_client

logger = logging.getLogger('signor')


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
}


_mechanism_map = {
    'catalytic activity': None, # Conversion
    'oxidoreductase activity': None,
    'transcriptional activation': None, # IncreaseAmount by tscript
    'transcriptional repression': None, # DecreaseAmount by tscript
    'Farnesylation': Farnesylation,
    'gtpase-activating protein': Gap,
    'deacetylation': Deacetylation,
    'demethylation': Demethylation,
    'dephosphorylation': Dephosphorylation,
    'destabilization': DecreaseAmount,
    'guanine nucleotide exchange factor': Gef,
    'acetylation': Acetylation,
    'binding': Complex,
    'cleavage': None, # Important!
    'desumoylation': Desumoylation,
    'deubiquitination': Deubiquitination,
    'glycosylation': Glycosylation,
    'hydroxylation': Hydroxylation,
    'neddylation': None, # Neddylation,
    'chemical activation': Activation,
    'chemical inhibition': Inhibition,
    'trimethylation': Methylation,
    'ubiquitination': Ubiquitination,
    'post transcriptional regulation': None,
    'relocalization': None, # Translocation,
    'small molecule catalysis': None,
    's-nitrosylation': None,
    'transcriptional regulation': None, # Need to know if up or down
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
    'down-regulates': Inhibition, # FIXME
    'down-regulates activity': Inhibition,
    'down-regulates quantity': DecreaseAmount,
    'down-regulates quantity by destabilization': DecreaseAmount,
    'down-regulates quantity by repression': DecreaseAmount,
    'form complex': Complex,
    'unknown': None,
    'up-regulates': Activation, # FIXME
    'up-regulates activity': Activation,
    'up-regulates quantity': IncreaseAmount,
    'up-regulates quantity by expression': IncreaseAmount,
    'up-regulates quantity by stabilization': IncreaseAmount
}


signor_default_path = join(dirname(__file__), '..', '..', 'data',
                          'all_data_23_09_17.csv')


class SignorProcessor(object):
    """Processor for Signor dataset, available at http://signor.uniroma2.it.

    See publication:

    Perfetto et al., "SIGNOR: a database of causal relationships between
    biological entities," Nucleic Acids Research, Volume 44, Issue D1, 4
    January 2016, Pages D548-D554. https://doi.org/10.1093/nar/gkv1048

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
    skipped_rows : list of SignorRow namedtuples
        List of rows where no mechanism statements were generated.
    skip_ctr : collections.Counter
        Counter listing the frequency of different MECHANISM types in the
        list of skipped rows.
    """
    def __init__(self, signor_csv=None, delimiter=';'):
        # Get generator over the CSV file
        if signor_csv:
            data_iter = read_unicode_csv(signor_csv, delimiter=';', skiprows=1)
            # Process into a list of SignorRow namedtuples
            # Strip off any funky \xa0 whitespace characters
            self._data = [SignorRow(*[f.strip() for f in r]) for r in data_iter]
        # If no CSV given, download directly from web
        else:
            res = requests.post('http://signor.uniroma2.it/download_entity.php',
                                data={'organism':'human', 'format':'csv',
                                      'submit':'Download'})
            if res.status_code == 200:
                csv_reader = csv.reader(StringIO(res.text), delimiter=delimiter,
                                        quoting=csv.QUOTE_MINIMAL)
                next(csv_reader) # Skip the header row
                self._data = [SignorRow(*[f.strip() for f in r])
                              for r in csv_reader]
            else:
                raise Exception('Could not download Signor data.')
        # Process into statements
        self.statements = []
        self.skipped_rows = []
        for row in self._data:
            row_stmts = self._process_row(row)
            self.statements.extend(row_stmts)
        # Skipped statements by type
        skip_ctr = Counter([row.MECHANISM for row in self.skipped_rows])
        self.skip_ctr = sorted([(k, v) for k, v in skip_ctr.items()],
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
            assert database in ('PUBCHEM', 'SIGNOR')
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
        return Evidence(source_api='SIGNOR', source_id=row.SIGNOR_ID,
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
            return stmts

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
                    stmts.append(ActiveForm(af_agent, 'activity', True))
                elif row.EFFECT.startswith('down'):
                    stmts.append(ActiveForm(af_agent, 'activity', False))
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
                    stmts.append(ActiveForm(af_agent, 'activity', True))
                elif row.EFFECT.startswith('down'):
                    stmts.append(ActiveForm(af_agent, 'activity', False))
        # For Complex statements, we create an ActiveForm with a BoundCondition.
        elif mech_stmt_type == Complex:
            # Complex
            stmts.append(mech_stmt_type([agent_a, agent_b], evidence=evidence))
            # ActiveForm
            af_agent = deepcopy(agent_b)
            af_bc_agent = deepcopy(agent_a)
            af_agent.bound_conditions = [BoundCondition(af_bc_agent, True)]
            if row.EFFECT.startswith('up'):
                stmts.append(ActiveForm(af_agent, 'activity', True))
            elif row.EFFECT.startswith('down'):
                stmts.append(ActiveForm(af_agent, 'activity', False))
        # Other mechanism statement types
        elif mech_stmt_type:
            stmts.append(mech_stmt_type(agent_a, agent_b, evidence=evidence))

        return stmts


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

"""
Known issues:
* The generic "up-regulates" effect type should be mapped to a generic up
  regulation rather than Activation/Inhibition, as it is currently.
* Mappings for SIGNOR families, complexes, etc.
* Whether Gef/Gap should produce additional statements to Act/Inh.
* ActiveForms representing effects on amounts (StateEffect)
"""
