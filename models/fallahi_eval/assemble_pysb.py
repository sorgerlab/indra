import re
import numpy
import pickle
from pysb.integrate import Solver
from indra.statements import *
from indra.mechlinker import MechLinker
import indra.tools.assemble_corpus as ac
from indra.databases import context_client, cbio_client, hgnc_client, \
                            uniprot_client
from indra.assemblers import PysbAssembler, IndexCardAssembler
from util import prefixed_pkl, pklload
from process_data import antibody_map, cell_lines, read_ccle_variants, \
                         drug_targets, drug_grounding, agent_from_gene_name

def assemble_pysb(stmts, data_genes, contextualize=False):
    # Filter the INDRA Statements to be put into the model
    stmts = ac.filter_by_type(stmts, Complex, invert=True)
    stmts = ac.filter_direct(stmts)
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    # Strip the extraneous supports/supported by here
    strip_supports(stmts)
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
    stmts = normalize_active_forms(ml.statements)
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
    # Save the Statements here
    ac.dump_statements(stmts, prefixed_pkl('pysb_stmts'))


    # Add drug target Statements
    drug_target_stmts = get_drug_target_statements()
    stmts += drug_target_stmts

    # Just generate the generic model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    with open(prefixed_pkl('pysb_model'), 'wb') as f:
        pickle.dump(model, f)

    # Run this extra part only if contextualize is set to True
    if not contextualize:
        return

    cell_lines_no_data = ['COLO858', 'K2', 'MMACSF', 'MZ7MEL', 'WM1552C']
    for cell_line in cell_lines:
        if cell_line not in cell_lines_no_data:
            stmtsc = contextualize_stmts(stmts, cell_line, data_genes)
        else:
            stmtsc = stmts
        pa = PysbAssembler()
        pa.add_statements(stmtsc)
        model = pa.make_model()
        if cell_line not in cell_lines_no_data:
            contextualize_model(model, cell_line, data_genes)
        ac.dump_statements(stmtsc, prefixed_pkl('pysb_stmts_%s' % cell_line))
        with open(prefixed_pkl('pysb_model_%s' % cell_line), 'wb') as f:
            pickle.dump(model, f)


def strip_supports(stmts):
    for stmt in stmts:
        stmt.supports = []
        stmt.supported_by = []


def get_drug_target_statements():
    stmts = []
    for drug, targets in drug_targets.items():
        for target in targets:
            target_agent = agent_from_gene_name(target)
            drug_agent = Agent(drug, db_refs=drug_grounding[drug])
            st = DecreaseAmount(drug_agent, target_agent)
            stmts.append(st)
    return stmts


def normalize_active_forms(stmts):
    af_stmts = ac.filter_by_type(stmts, ActiveForm)
    relevant_af_stmts = []
    for stmt in af_stmts:
        if (not stmt.agent.mods) and (not stmt.agent.mutations):
            continue
        relevant_af_stmts.append(stmt)
    print('%d relevant ActiveForms' % len(relevant_af_stmts))
    non_af_stmts = ac.filter_by_type(stmts, ActiveForm, invert=True)
    af_stmts = ac.run_preassembly(relevant_af_stmts)
    stmts = af_stmts + non_af_stmts
    return stmts


def contextualize_stmts(stmts, cell_line, genes):
    """Contextualize model at the level of INDRA Statements."""
    to_remove = []
    cell_line_ccle = cell_line + '_SKIN'
    # Remove genes with CNA = -2
    print('Contextualize by CNA')
    cna = cbio_client.get_ccle_cna(genes, [cell_line_ccle])[cell_line_ccle]
    for gene in genes:
        if cna.get(gene) == -2:
            to_remove.append(gene)
            print('To remove CNA: %s' % gene)
    # Remove genes with transcripts in bottom 5%
    print('Contextualize by mRNA')
    mrna = cbio_client.get_ccle_mrna(genes, [cell_line_ccle])[cell_line_ccle]
    mrna_vals = [v for v in mrna.values() if v]
    if mrna_vals:
        thresh = numpy.percentile(mrna_vals, 5.0)
        for gene, val in mrna.items():
            if val and val < thresh:
                to_remove.append(gene)
                print('To remove mRNA: %s' % gene)
    # Remove genes with nonsense mutations
    print('Contextualize by nonsense mutations')
    variants = read_ccle_variants(genes)
    to_remove_nonsense = list(variants['nonsense'][cell_line_ccle].keys())
    if to_remove_nonsense:
        print('To remove nonsense: %s' % ', '.join(to_remove_nonsense))
        to_remove += to_remove_nonsense
    # Remove Statements for these genes
    new_stmts = []
    for stmt in stmts:
        any_to_remove = False
        for agent in stmt.agent_list():
            if agent is not None and agent.name in to_remove:
                any_to_remove = True
                break
        if not any_to_remove:
            new_stmts.append(stmt)

    # Remove Statements with irrelevant mutations
    print('Contextualize by missense mutations')
    mutations = variants['missense'][cell_line_ccle]
    muts_to_use = {}
    for gene, mut_list in mutations.items():
        muts_to_use[gene] = []
        for mut in mut_list:
            muts_to_use[gene].append(mut)
    stmts = ac.filter_mutation_status(new_stmts, mutations, [])
    return stmts


def contextualize_model(model, cell_line, genes):
    """Contextualize model at the level of a PySB model."""
    # Here we just make a PysbAssembler to be able
    # to apply set_context on the model being passed in
    model.name = cell_line
    cell_line_ccle = cell_line + '_SKIN'
    pa = PysbAssembler()
    pa.model = model
    pa.set_context(cell_line_ccle)

    # Set initial conditions for missense mutations
    variants = read_ccle_variants(genes)
    mutations = variants['missense'][cell_line_ccle]
    for gene, mut_list in mutations.items():
        for fres, loc, tres in mut_list:
            site_name = fres + loc
            for ic in model.initial_conditions:
                if ic[0].monomer_patterns[0].monomer.name == gene:
                    sc = ic[0].monomer_patterns[0].site_conditions
                    if site_name in sc:
                        sc[site_name] = tres

    return pa.model


def get_mod_whitelist():
    """Return dict of modifications that are relevant for the data."""
    mod_whitelist = {}
    for members in antibody_map.values():
        for gene_name, phos_sites in members.items():
            for residue, position in phos_sites:
                if gene_name not in mod_whitelist:
                    mod_whitelist[gene_name] = []
                entry = ('phosphorylation', residue, position)
                mod_whitelist[gene_name].append(entry)
    return mod_whitelist


def print_initial_conditions(cell_lines, gene_names, fname):
    """Print a table with initial states/conditions for each cell line."""
    fh = open(fname, 'wt')
    genes_per_model = {}
    # First we collect the state and initial value of proteins
    # in each cell line into a dict
    for cell_line in cell_lines:
        genes_per_model[cell_line] = {}
        model = pklload('pysb_model_%s' % cell_line)
        inits = sorted(model.initial_conditions,
            key=lambda x: x[0].monomer_patterns[0].monomer.name)
        # A list of "basic" states, deviation from which indicates
        # a "special" state like a mutated form
        basic_states = ['n', 'u', 'WT', 'inactive', None]
        # Iterate over all the initial conditions to collect them
        # in the genes_per_model datastructure
        for cplx, initial in inits:
            name = cplx.monomer_patterns[0].monomer.name
            extra_states = []
            for site, state in cplx.monomer_patterns[0].site_conditions.items():
                if site == 'loc':
                    continue
                if state not in basic_states:
                    extra_states.append('%s=%s' % (site, state))

            val_str = '%d' % initial.value
            if extra_states:
                val_str += ' (%s)' % ', '.join(extra_states)
            genes_per_model[cell_line][name] = val_str
    # Print the file while skipping proteins that aren't in any of
    # the cell line specific models.
    header = 'Protein\t' + '\t'.join(cell_lines) + '\n'
    fh.write(header)
    for gene in gene_names:
        line_vals = [gene]
        no_val = True
        for cell_line in cell_lines:
            val = genes_per_model[cell_line].get(gene)
            if not val:
                val = '-'
            else:
                no_val = False
            line_vals.append(val)
        line = '\t'.join(line_vals)
        if not no_val:
            fh.write(line + '\n')
    fh.close()


def process_kappa_dead_rules(kasa_output):
    dead_rules = []
    for line in kasa_output.split('\n'):
        match = re.match('rule ([\d]+): ([a-zA-Z0-9_]+) will never be applied.',
                         line.strip())
        if match:
            dead_rules.append(match.groups()[1])
    return dead_rules


def remove_kappa_dead_rules(stmts, model, dead_rules):
    # FIXME: we should probably check that a statement we remove has all its
    # generated rules recognized as dead. If it has at least one live rule
    # coming from it, we shouldn't remove it. But the dead rules should still
    # be removed somehow from the final model.
    dead_uuids  = set()
    for rule in dead_rules:
        for ann in model.annotations:
            if ann.subject == rule and ann.predicate == 'from_indra_statement':
                dead_uuids.add(ann.object)
    all_uuids = {stmt.uuid for stmt in stmts}
    live_uuids = all_uuids - dead_uuids
    stmts = ac.filter_uuid_list(stmts, live_uuids)
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    return stmts, model
