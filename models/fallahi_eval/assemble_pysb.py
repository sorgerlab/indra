import pickle
from indra.statements import *
from indra.mechlinker import MechLinker
import indra.tools.assemble_corpus as ac
from indra.assemblers import PysbAssembler, IndexCardAssembler
from process_data import antibody_map, cell_lines
from indra.databases import context_client
from util import prefixed_pkl

def assemble_pysb(stmts, data_genes, contextualize=False):
    # Filter the INDRA Statements to be put into the model
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
    ml.gather_modifications()
    ml.reduce_modifications()
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

    # Just generate the generic model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    ac.dump_statements(stmts, prefixed_pkl('pysb_stmts'))
    with open(prefixed_pkl('pysb_model'), 'wb') as f:
        pickle.dump(model, f)

    # Run this extra part only if contextualize is set to True
    if not contextualize:
        return

    cell_lines_no_data = ['COLO858', 'K2', 'MMACSF', 'MZ7MEL', 'WM1552C']
    for cell_line in cell_lines:
        if cell_line not in cell_lines_no_data:
            stmtsc = contextualize_stmts(stmts, 'C32', data_genes)
        else:
            stmtsc = stmts
        pa = PysbAssembler()
        pa.add_statements(stmtsc)
        model = pa.make_model()
        if cell_line not in cell_lines_no_data:
            contextualize_model(model, cell_line)
        ac.dump_statements(stmtsc, prefixed_pkl('pysb_stmts_%s' % cell_line))
        with open(prefixed_pkl('pysb_model_%s' % cell_line), 'wb') as f:
            pickle.dump(model, f)


def contextualize_stmts(stmts, cell_line, genes):
    mutations = context_client.get_mutations(genes, [cell_line])
    mutations = mutations.get(cell_line)
    muts_to_use = {}
    for gene, mut_list in mutations.items():
        muts_to_use[gene] = []
        for mut in mut_list:
            try:
                from_aa = mut[0]
                to_aa = mut[-1]
                pos = mut[1:-1]
                muts_to_use.append((from_aa, pos, to_aa))
            except Exception:
                print(mut)
    stmts = ac.filter_mutation_status(stmts, mutations, [])
    return stmts


def contextualize_model(model, cell_line):
    # Here we just make a PysbAssembler to be able
    # to apply set_context on the model being passed in
    if not cell_line.endswith('_SKIN'):
        cell_line = cell_line + '_SKIN'
    pa = PysbAssembler()
    pa.model = model
    pa.set_context(cell_line)
    return pa.model

def get_mod_whitelist():
    mod_whitelist = {}
    for members in antibody_map.values():
        for gene_name, phos_sites in members.items():
            for residue, position in phos_sites:
                if gene_name not in mod_whitelist:
                    mod_whitelist[gene_name] = []
                entry = ('phosphorylation', residue, position)
                mod_whitelist[gene_name].append(entry)
    return mod_whitelist
