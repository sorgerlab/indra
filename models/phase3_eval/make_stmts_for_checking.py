import pickle
import itertools
from collections import defaultdict
import process_data as pd
from indra.databases import hgnc_client
from read_phosphosite import read_phosphosite
from indra.statements import Agent, Dephosphorylation, Phosphorylation
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies

drug_col = 'Sample Description (drug abbre. | dose or time-point)'

def get_target_agent(target):
    target_hgnc_id = hgnc_client.get_hgnc_id(target)
    target_up_id = hgnc_client.get_uniprot_id(target_hgnc_id)
    target_agent = Agent(target, db_refs={'HGNC': target_hgnc_id,
                                          'UP': target_up_id})
    return target_agent

def get_drug_data(data, tx):
    drug_tx_data = data['protein'][data['protein'][drug_col] == tx]
    return drug_tx_data

def preassemble_stmts(stmts):
    pa_dict = {} # Preassembled dict (convert to regular dict because
                 # defaultdict with lambda can't be preassembled
    for drug_name, ab_dict in stmts.items():
        pa_ab_dict = {}
        for ab_name, stmt_list in ab_dict.items():
            pa = Preassembler(hierarchies)
            pa.add_statements(stmt_list)
            pa.combine_duplicates()
            pa_ab_dict[ab_name] = pa.unique_stmts
        pa_dict[drug_name] = pa_ab_dict
    return pa_dict

def get_phospho_stmts(target_agent, phosphoforms, fold_change):
    stmts = []
    # Each entry in this list is an Agent with db_refs filled in
    # and associated mod conditions for the phosphorylated sites
    for psf in phosphoforms:
        psf_agent = Agent(psf.name, db_refs=psf.db_refs)
        for mod in psf.mods:
            # Create a Phosphorylation statement corresponding to
            # this drug/Ab pair
            if fold_change < 1:
                stmt = Phosphorylation(target_agent, psf_agent,
                                       mod.residue, mod.position)
            else:
                stmt = Dephosphorylation(target_agent, psf_agent,
                                         mod.residue, mod.position)
            stmts.append(stmt)
    return stmts

def make_stmts(data, ab_agents, drug_ab_combs=None, thresh=None):
    if drug_ab_combs is None:
        drug_tx = pd.get_single_drug_treatments(data)
        antibodies = pd.get_phos_antibodies(data)
        drug_ab_combs = itertools.product(drug_tx, antibodies)

    dec_thresh, inc_thresh = thresh if thresh is not None else (1, 1)
    drug_targets = pd.get_drug_targets()
    stmts = defaultdict(lambda: defaultdict(list))
    values = defaultdict(dict)

    for tx, ab in drug_ab_combs:
        # Get the drug name
        drug_name = tx.split('|')[0]
        targets = drug_targets[drug_name]
        drug_tx_data = get_drug_data(data, tx)
        for target in targets:
            # Create an agent for the drug target
            target_agent = get_target_agent(target)
            fold_change = drug_tx_data[ab].values[0]
            if fold_change < dec_thresh or fold_change > inc_thresh:
                phosphoforms = ab_agents[ab]
                phos_stmts = get_phospho_stmts(target_agent, phosphoforms,
                                               fold_change)
                stmts[drug_name][ab] += phos_stmts
                values[drug_name][ab] = fold_change
    return stmts, values

def get_eval_drug_ab_combs(data):
    drug_tx = data['prediction'][drug_col]
    antibodies = data['prediction']['Antibody ID']
    drug_ab_combs = zip(drug_tx, antibodies)
    return drug_ab_combs

def run(dec_thresh=0.8, inc_thresh=1.2):
    data = pd.read_data(pd.data_file)
    ab_agents = read_phosphosite('sources/annotated_kinases_v5.csv')[1]

    # If filtering is to be done based on thresholds only,
    # set this to None
    drug_ab_combs = get_eval_drug_ab_combs(data)

    stmts, values = make_stmts(data, ab_agents, drug_ab_combs=drug_ab_combs,
                               thresh=[dec_thresh, inc_thresh])

    # Now, preassemble the statements to remove duplicates
    pa_dict = preassemble_stmts(stmts)

    with open('data_stmts.pkl', 'wb') as f:
        pickle.dump((pa_dict, values), f, protocol=2)

    return (stmts, values)

if __name__ == '__main__':
    run()
