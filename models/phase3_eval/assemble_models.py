import pickle
import itertools
from copy import deepcopy
from os.path import join as pjoin
from pysb.export.kappa import KappaExporter
from indra.tools import assemble_corpus as ac
from indra.tools.gene_network import GeneNetwork
from indra.assemblers import IndexCardAssembler
from indra.statements import *
from indra.assemblers import PysbAssembler, SifAssembler, CxAssembler

import process_data
from read_phosphosite import read_phosphosite

def build_prior(genes, out_file):
    gn = GeneNetwork(genes, 'korkut')
    stmts = gn.get_statements(filter=False)
    ac.dump_statements(stmts, out_file)
    return stmts

def assemble_index_cards(stmts):
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
            ica.save_model('output/index_cards/index_card_%d.json' % counter)
            counter += 1

def assemble_pysb(stmts, data_genes, out_file):
    """Return an assembled PySB model."""
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_gene_list(stmts, data_genes, 'all')
    stmts = ac.reduce_activities(stmts)
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()
    pa.set_context('SKMEL28_SKIN')
    pa.save_model(out_file)
    return model

def assemble_sif(stmts, data, out_file):
    """Return an assembled SIF."""
    # Filter for high-belief statements
    stmts = ac.filter_belief(stmts, 0.99)
    # Filter for Activation / Inhibition
    stmts_act = ac.filter_by_type(stmts, Activation)
    stmts_inact = ac.filter_by_type(stmts, Inhibition)
    stmts = stmts_act + stmts_inact
    # Get the drugs inhibiting their targets as INDRA
    # statements
    def get_drug_statements():
        drug_targets = process_data.get_drug_targets()
        drug_stmts = []
        for dn, tns in drug_targets.items():
            da = Agent(dn + ':Drugs')
            for tn in tns:
                ta = Agent(tn)
                drug_stmt = Inhibition(da, ta)
                drug_stmts.append(drug_stmt)
        return drug_stmts
    drug_stmts = get_drug_statements()
    stmts = stmts + drug_stmts
    # Because of a bug in CNO, node names containing AND
    # need to be replaced
    def rename_and_nodes(st):
        for s in st:
            for a in s.agent_list():
                if a is not None:
                    if a.name.find('AND') != -1:
                        a.name = a.name.replace('AND', 'A_ND')
    rename_and_nodes(stmts)
    # Rewrite statements to replace genes with their corresponding
    # antibodies when possible
    stmts = rewrite_ab_stmts(stmts, data)
    def filter_ab_edges(st, policy='all'):
        st_out = []
        for s in st:
            if policy == 'all':
                all_ab = True
                for a in s.agent_list():
                    if a is not None:
                        if a.name.find('_p') == -1 and \
                            a.name.find('Drugs') == -1:
                            all_ab = False
                            break
                if all_ab:
                    st_out.append(s)
            elif policy == 'one':
                any_ab = False
                for a in s.agent_list():
                    if a is not None and a.name.find('_p') != -1:
                        any_ab = True
                        break
                if any_ab:
                    st_out.append(s)
        return st_out
    stmts = filter_ab_edges(stmts, 'one')
    print(len(stmts))
    # Make the SIF model
    sa = SifAssembler(stmts)
    sa.make_model(use_name_as_key=True)
    sif_str = sa.print_model()
    with open(out_file, 'wb') as fh:
        fh.write(sif_str.encode('utf-8'))
    # Make the MIDAS data file used for training the model
    midas_data = process_data.get_midas_data(data)
    return sif_str

def assemble_cx(stmts, out_file):
    """Return a CX assembler."""
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.strip_agent_context(stmts)
    ca = CxAssembler()
    ca.add_statements(stmts)
    model = ca.make_model()
    ca.save_model(out_file)
    return ca

def get_prior_genes(fname):
    """Get the list of prior genes."""
    with open(fname, 'rt') as fh:
        genes = fh.read().strip().split('\n')
        return genes

def rewrite_ab_stmts(stmts_in, data):
    """Replace corresponding Agent names to AB names.

    This is used for logical model building.
    """
    # Build AB <-> Gene maps
    ab_phos = process_data.get_phos_antibodies(data)
    ab_dict = {}
    ab_dict_rev = {}
    for ab in ab_phos:
        up_ids = process_data.get_antibody_genes(data, ab)
        for up_id in up_ids:
            agent = process_data.get_agent_from_upid(up_id)
            if ab in ab_dict:
                ab_dict[ab].append(agent)
            else:
                ab_dict[ab] = [agent]
            if agent.name in ab_dict_rev:
                ab_dict_rev[agent.name].append(ab)
            else:
                ab_dict_rev[agent.name] = [ab]
    stmts_out = []
    for stmt in stmts_in:
        new_agents = []
        for agent in stmt.agent_list():
            if agent is not None:
                ab_names = ab_dict_rev.get(agent.name)
                if not ab_names:
                    new_agents.append([agent])
                else:
                    new_agents.append([Agent(ab) for ab in ab_names])
            else:
                new_agents.append([None])
        for na in itertools.product(*new_agents):
            stmt_new = deepcopy(stmt)
            stmt_new.set_agent_list(na)
            stmts_out.append(stmt_new)
    return stmts_out

if __name__ == '__main__':
    outf = 'output/'
    data = process_data.read_data(process_data.data_file)
    data_genes = process_data.get_all_gene_names(data)
    reassemble = False
    if not reassemble:
        stmts = ac.load_statements(pjoin(outf, 'top_level.pkl'))
    else:
        #prior_stmts = build_prior(data_genes, pjoin(outf, 'prior.pkl'))
        prior_stmts = ac.load_statements(pjoin(outf, 'prior.pkl'))
        prior_stmts = ac.map_grounding(prior_stmts, dump_pkl=pjoin(outf, 'gmapped_prior.pkl'))
        reading_stmts = ac.load_statements(pjoin(outf, 'phase3_stmts.pkl'))
        reading_stmts = ac.map_grounding(reading_stmts, dump_pkl=pjoin(outf, 'gmapped_reading.pkl'))
        stmts = prior_stmts + reading_stmts

        stmts = ac.filter_grounded_only(stmts)
        stmts = ac.filter_genes_only(stmts, specific_only=False)
        stmts = ac.filter_human_only(stmts)
        stmts = ac.expand_families(stmts)
        stmts = ac.filter_gene_list(stmts, data_genes, 'one')
        stmts = ac.map_sequence(stmts, dump_pkl=pjoin(outf, 'smapped.pkl'))
        stmts = ac.run_preassembly(stmts, dump_pkl=pjoin(outf, 'top_level.pkl'))

    assemble_models = []
    assemble_models.append('sif')
    assemble_models.append('pysb')
    assemble_models.append('cx')

    ### PySB assembly
    if 'pysb' in assemble_models:
        pysb_model = assemble_pysb(stmts, data_genes,
                                   pjoin(outf, 'korkut_model_pysb.py'))
        ke = KappaExporter(pysb_model)
        with open(pjoin(outf, 'korkut_model.ka'), 'wb') as fh:
            fh.write(ke.export().encode('utf-8'))
    ### SIF assembly
    if 'sif' in assemble_models:
        sif_str = assemble_sif(stmts, data, pjoin(outf, 'PKN-korkut_one_ab.sif'))
    ### CX assembly
    if 'cx' in assemble_models:
        cxa = assemble_cx(stmts, pjoin(outf, 'korkut_full_high_belief.cx'))
