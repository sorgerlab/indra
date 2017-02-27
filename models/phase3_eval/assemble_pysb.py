from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from os.path import join as pjoin
import os.path
from pysb import Observable, bng
from pysb.export.kappa import KappaExporter
from indra.util import read_unicode_csv
from indra.assemblers import PysbAssembler, IndexCardAssembler
from indra.mechlinker import MechLinker
from indra.statements import *
import indra.tools.assemble_corpus as ac
from read_phosphosite import read_phosphosite

def assemble_pysb(stmts, data_genes, out_file):
    """Return an assembled PySB model."""
    base_file, _ = os.path.splitext(out_file)
    # IF YOU DON'T WANT OT RERUN THE PREPROCESSING, LOAD FROM
    # A PICKLE HERE
    # stmts = ac.load_statements('pysb_assembler_input.pkl')
    stmts = preprocess_stmts(stmts, data_genes)

    # This is the "final" set of statements going into the assembler so it
    # makes sense to cache these.
    # This is also the point where index cards can be generated
    ac.dump_statements(stmts, '%s.pkl' % base_file)

    # Assemble model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model()
    # Set context
    set_context(pa)
    # Add observables
    add_observables(pa.model)
    pa.save_model(out_file)
    pa.export_model('kappa', '%s.ka' % base_file)
    return pa.model

def generate_equations(model, pkl_cache):
    bng.generate_equations(model, verbose=True)
    with open(pkl_cache, 'w') as fh:
        pickle.dump(model, fh)

def set_context(pa):
    pa.set_context('SKMEL28_SKIN')
    # Set BRAF V600E
    for ic in pa.model.initial_conditions:
        if str(ic[0]).startswith('BRAF'):
            ic[0].monomer_patterns[0].site_conditions['V600'] = 'E'

def add_observables(model):
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
    o = Observable(b'EIF4EBP1p', model.monomers['EIF4EBP1'](phospho='p'))
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
    o = Observable(b'RPS6KA1p', model.monomers['RPS6KA1'](phospho='p'))
    model.add_component(o)
    o = Observable(b'RPS6KB1p', model.monomers['RPS6KB1'](T359='p'))
    model.add_component(o)
    o = Observable(b'PDPK1p', model.monomers['PDPK1'](S241='p'))
    model.add_component(o)
    o = Observable(b'PTK2p', model.monomers['PTK2'](Y397='p'))
    model.add_component(o)
    o = Observable(b'STAT3p', model.monomers['STAT3'](S727='p'))
    model.add_component(o)
    o = Observable(b'IRS1p', model.monomers['IRS1'](S307='p'))
    model.add_component(o)
    o = Observable(b'ESR1p', model.monomers['ESR1p'](S118='p'))
    model.add_component(o)


def get_mod_whitelist():
    mod_whitelist = {}
    _, ab_map = read_phosphosite('annotated_kinases_v4.csv')
    for k, v in ab_map.items():
        for agent in v:
            res_pos = (agent.mods[0].residue, agent.mods[0].position)
            try:
                mod_whitelist[agent.name].append(res_pos)
            except KeyError:
                mod_whitelist[agent.name] = [res_pos]
    return mod_whitelist


def preprocess_stmts(stmts, data_genes):
    # Filter the INDRA Statements to be put into the model
    stmts = filter_mutation_status(stmts, {'BRAF': [('V', '600', 'E')]}, ['PTEN'])
    stmts = ac.filter_by_type(stmts, Complex, invert=True)
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_gene_list(stmts, data_genes, 'all')
    stmts = filter_enzyme_kinase(stmts)
    stmts = filter_mod_nokinase(stmts)
    stmts = filter_transcription_factor(stmts)
    # Simplify activity types
    ml = MechLinker(stmts)
    ml.get_explicit_activities()
    ml.reduce_activities()
    af_stmts = ac.filter_by_type(ml.statements, ActiveForm)
    non_af_stmts = ac.filter_by_type(ml.statements, ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(af_stmts)
    stmts = af_stmts + non_af_stmts
    # Replace activations when possible
    ml = MechLinker(stmts)
    ml.get_explicit_activities()
    ml.replace_activations()
    # Require active forms
    ml.require_active_form()
    # Remove inconsequential PTMs
    ml.statements = filter_inconsequential_ptms(ml.statements,
                                                get_mod_whitelist())
    stmts = ml.statements
    return stmts


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

# These filters are currently in a very basic implementation stage
# where they work well on the statements encountered in the scope of this
# evaluation. They could be generalized and refactored to become
# part of the core of INDRA.
# Note that implicitly these filters also filter out statements in which
# the subject is None.
def filter_enzyme_kinase(stmts_in):
    kinase_table = read_unicode_csv('../../indra/resources/kinases.tsv',
                                    delimiter='\t')
    gene_names = [lin[1] for lin in list(kinase_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, Phosphorylation):
            if st.enz is not None:
                if st.enz.name in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    return stmts_out

def filter_mod_nokinase(stmts_in):
    kinase_table = read_unicode_csv('../../indra/resources/kinases.tsv',
                                    delimiter='\t')
    gene_names = [lin[1] for lin in list(kinase_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, Modification) and not \
           isinstance(st, Phosphorylation):
            if st.enz is not None:
                if st.enz.name not in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    return stmts_out

def filter_transcription_factor(stmts_in):
    tf_table = \
        read_unicode_csv('../../indra/resources/transcription_factors.csv')
    gene_names = [lin[1] for lin in list(tf_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, RegulateAmount):
            if st.subj is not None:
                if st.subj.name in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    return stmts_out

def filter_inconsequential_ptms(stmts_in, whitelist):
    states_used = whitelist
    for stmt in stmts_in:
        for agent in stmt.agent_list():
            if agent is not None:
                if agent.mods:
                    for mc in agent.mods:
                        res_pos = (mc.residue, mc.position)
                        try:
                            states_used[agent.name].append(res_pos)
                        except KeyError:
                            states_used[agent.name] = [res_pos]
    stmts_out = []
    for stmt in stmts_in:
        skip = False
        if isinstance(stmt, Modification):
            if stmt.residue is not None or stmt.position is not None:
                res_pos = (stmt.residue, stmt.position)
                used = states_used.get(stmt.sub.name, [])
                if res_pos not in used:
                    skip = True
        if not skip:
            stmts_out.append(stmt)
        #else:
        #    print('Removed: %s' % stmt)
    print('Total statements remaining: %d' % len(stmts_out))
    return stmts_out

def filter_mutation_status(stmts_in, mutations, deletions):
    stmts_out = []
    for stmt in stmts_in:
        skip = False
        for agent in stmt.agent_list():
            if agent is not None and agent.name in deletions:
                skip = True
            if agent is not None and agent.mutations:
                muts = mutations.get(agent.name, [])
                for mut in agent.mutations:
                    mut_tup = (mut.residue_from, mut.position, mut.residue_to)
                    if mut_tup not in muts:
                        skip = True
            if skip:
                break
        if not skip:
            stmts_out.append(stmt)
        #else:
        #    print('Removed: %s' % stmt)
    print('Total statements remaining: %d' % len(stmts_out))
    return stmts_out
