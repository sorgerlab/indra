from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import pickle
import os.path
from os.path import join as pjoin
from pysb import Observable, ReactionPattern, ComplexPattern, bng
from pysb.export.kappa import KappaExporter
from indra.statements import *
from indra.mechlinker import MechLinker
import indra.tools.assemble_corpus as ac
from indra.assemblers import PysbAssembler, IndexCardAssembler
from read_phosphosite import read_phosphosite

def assemble_pysb(stmts, data_genes, out_file):
    """Return an assembled PySB model."""
    base_file, _ = os.path.splitext(out_file)
    #stmts = ac.load_statements('%s.pkl' % base_file)
    stmts = preprocess_stmts(stmts, data_genes)

    # This is the "final" set of statements going into the assembler so it
    # makes sense to cache these.
    # This is also the point where index cards can be generated
    ac.dump_statements(stmts, '%s_before_pa.pkl' % base_file)
    assemble_index_cards(stmts, 'output/index_cards')

    # Assemble model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model(reverse_effects=True)
    #ac.dump_statements(pa.statements, '%s_after_pa.pkl' % base_file)
    # Set context
    set_context(pa)
    # Add observables
    add_observables(pa.model)
    pa.save_model(out_file)
    #pa.export_model('kappa', '%s.ka' % base_file)
    return pa.model

def preprocess_stmts(stmts, data_genes):
    # Filter the INDRA Statements to be put into the model
    stmts = ac.filter_mutation_status(stmts,
                                      {'BRAF': [('V', '600', 'E')]}, ['PTEN'])
    stmts = ac.filter_by_type(stmts, Complex, invert=True)
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_gene_list(stmts, data_genes, 'all')
    stmts = ac.filter_enzyme_kinase(stmts)
    stmts = ac.filter_mod_nokinase(stmts)
    stmts = ac.filter_transcription_factor(stmts)
    # Simplify activity types
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.reduce_activities()
    af_stmts = ac.filter_by_type(ml.statements, ActiveForm)
    non_af_stmts = ac.filter_by_type(ml.statements, ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(af_stmts)
    stmts = af_stmts + non_af_stmts
    # Replace activations when possible
    ml = MechLinker(stmts)
    ml.gather_explicit_activities()
    ml.replace_activations()
    # Require active forms
    ml.require_active_forms()
    num_stmts = len(ml.statements)
    while True:
        # Remove inconsequential PTMs
        ml.statements = ac.filter_inconsequential_mods(ml.statements,
                                                       get_mod_whitelist())
        ml.statements = ac.filter_inconsequential_acts(ml.statements,
                                                       get_mod_whitelist())
        if num_stmts <= len(ml.statements):
            break
        num_stmts = len(ml.statements)
    stmts = ml.statements
    return stmts


def set_context(pa):
    pa.set_context('SKMEL28_SKIN')
    # Set BRAF V600E
    for ic in pa.model.initial_conditions:
        if str(ic[0]).startswith('BRAF'):
            ic[0].monomer_patterns[0].site_conditions['V600'] = 'E'


def generate_equations(model, pkl_cache):
    bng.generate_equations(model, verbose=True)
    with open(pkl_cache, 'w') as fh:
        pickle.dump(model, fh)


def add_observables(model):
    _, ab_map = read_phosphosite('sources/annotated_kinases_v5.csv')
    for ab_name, agents in ab_map.items():
        patterns = []
        for agent in agents:
            try:
                monomer = model.monomers[agent.name]
            except KeyError:
                continue
            mc = agent.mods[0]
            site_names = ['phospho', mc.residue]
            if mc.position is not None:
                site_names.append(mc.residue + mc.position)
            for site_name in site_names:
                try:
                    pattern = monomer(**{site_name: 'p'})
                    patterns.append(ComplexPattern([pattern], None))
                except Exception:
                    pass
        if patterns:
            obs_name = ab_name
            if not re.match(r'[_a-z][_a-z0-9]*\Z', obs_name, re.IGNORECASE):
                obs_name = obs_name.replace('-', '_')
            if not re.match(r'[_a-z][_a-z0-9]*\Z', obs_name, re.IGNORECASE):
                obs_name = 'p' + obs_name
            obs_name = obs_name.encode('utf-8')
            o = Observable(obs_name, ReactionPattern(patterns))
            model.add_component(o)
    '''
    o = Observable(b'MAPK1p', model.monomers['MAPK1'](T185='p', Y187='p'))
    model.add_component(o)
    o = Observable(b'MAPK3p', model.monomers['MAPK3'](T202='p', Y204='p'))
    model.add_component(o)
    o = Observable(b'MAPK14p', model.monomers['MAPK14'](T180='p'))
    model.add_component(o)
    o = Observable(b'GSK3Ap', model.monomers['GSK3A'](S21='p'))
    model.add_component(o)
    o = Observable(b'GSK3Bp', model.monomers['GSK3B'](S9='p'))
    model.add_component(o)
    o = Observable(b'RPS6pS235', model.monomers['RPS6'](S235='p'))
    model.add_component(o)
    o = Observable(b'RPS6pS240', model.monomers['RPS6'](S240='p'))
    model.add_component(o)
    o = Observable(b'EIF4EBP1p', model.monomers['EIF4EBP1'](S65='p'))
    model.add_component(o)
    o = Observable(b'JUNp', model.monomers['JUN'](S73='p'))
    model.add_component(o)
    o = Observable(b'FOXO3p', model.monomers['FOXO3'](S315='p'))
    model.add_component(o)
    o = Observable(b'AKT1p', model.monomers['AKT1'](S473='p'))
    model.add_component(o)
    o = Observable(b'AKT2p', model.monomers['AKT2'](S474='p'))
    model.add_component(o)
    o = Observable(b'AKT3p', model.monomers['AKT3'](phospho='p'))
    model.add_component(o)
    o = Observable(b'ELK1p', model.monomers['ELK1'](S383='p'))
    model.add_component(o)
    o = Observable(b'RB1p', model.monomers['RB1'](S807='p'))
    model.add_component(o)
    o = Observable(b'RPS6KA1p', model.monomers['RPS6KA1'](T359='p'))
    model.add_component(o)
    o = Observable(b'RPS6KB1p', model.monomers['RPS6KB1'](phospho='p'))
    model.add_component(o)
    o = Observable(b'PDPK1p', model.monomers['PDPK1'](S241='p'))
    model.add_component(o)
    o = Observable(b'PTK2p', model.monomers['PTK2'](Y397='p'))
    model.add_component(o)
    o = Observable(b'STAT3p', model.monomers['STAT3'](S727='p'))
    model.add_component(o)
    o = Observable(b'IRS1p', model.monomers['IRS1'](S307='p'))
    model.add_component(o)
    o = Observable(b'ESR1p', model.monomers['ESR1'](S118='p'))
    model.add_component(o)
    '''


def get_mod_whitelist():
    mod_whitelist = {}
    _, ab_map = read_phosphosite('sources/annotated_kinases_v5.csv')
    for k, v in ab_map.items():
        for agent in v:
            mod = ('phosphorylation', agent.mods[0].residue,
                   agent.mods[0].position)
            try:
                mod_whitelist[agent.name].append(mod)
            except KeyError:
                mod_whitelist[agent.name] = [mod]
            # Add generic mods to make sure we keep them
            mod = ('phosphorylation', agent.mods[0].residue, None)
            mod_whitelist[agent.name].append(mod)
            mod = ('phosphorylation', None, None)
            mod_whitelist[agent.name].append(mod)
    return mod_whitelist


def assemble_index_cards(stmts, out_folder):
    counter = 1
    for st in stmts:
        if isinstance(st, Modification) and st.enz is None:
            continue
        pmids = [ev.pmid for ev in st.evidence if ev.pmid is not None]
        if pmids:
            pmids = ','.join(['PMID%s' % pm for pm in list(set(pmids))])
        else:
            pmids = 'N/A'
        ica = IndexCardAssembler([st], pmc_override=pmids)
        ica.make_model()
        if ica.cards:
            ica.save_model(pjoin(out_folder, 'index_card_%d.json' % counter))
            counter += 1
