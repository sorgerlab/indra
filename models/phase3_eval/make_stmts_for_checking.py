import pickle
import itertools
from collections import defaultdict
import process_data as pd
import numpy as np
from indra.databases import hgnc_client
from read_phosphosite import read_phosphosite
from indra.statements import Agent, Dephosphorylation, Phosphorylation, \
                             IncreaseAmount, DecreaseAmount
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

def get_observed_stmts(target_agent, observed_agent_forms, fold_change):
    stmts = []
    for obsf in observed_agent_forms:
        obs_agent = Agent(obsf.name, db_refs=obsf.db_refs)
        # If the agent has a modification then we make Modification
        # statements to check
        if obsf.mods:
            for mod in obsf.mods:
                # Create a Phosphorylation statement corresponding to
                # this drug/Ab pair
                if fold_change < 1:
                    stmt = Phosphorylation(target_agent, obs_agent,
                                           mod.residue, mod.position)
                else:
                    stmt = Dephosphorylation(target_agent, obs_agent,
                                             mod.residue, mod.position)
        # Otherwise the observed change is in protein amounts so we make
        # RegulateAmount statements
        else:
            if fold_change < 1:
                stmt = IncreaseAmount(target_agent, obs_agent)
            else:
                stmt = DecreaseAmount(target_agent, obs_agent)
        stmts.append(stmt)
    return stmts

def make_stmts(data, ab_agents, drug_ab_combs=None, thresh=None):
    if drug_ab_combs is None:
        drug_tx = pd.get_single_drug_treatments(data)
        antibodies = pd.get_all_antibodies(data)
        drug_ab_combs = itertools.product(drug_tx, antibodies)

    dec_thresh, inc_thresh = np.log2(thresh if thresh is not None else (1, 1))
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
            fold_change = np.log2(drug_tx_data[ab].values[0])
            if fold_change < dec_thresh or fold_change > inc_thresh:
                observed_agent_forms = ab_agents[ab]
                obs_stmts = get_observed_stmts(target_agent,
                                               observed_agent_forms,
                                               fold_change)
                stmts[drug_name][ab] += obs_stmts
                values[drug_name] = {k: np.log2(list(v.values())[0])
                                     for k, v in
                                     drug_tx_data.iloc[:,2:].to_dict().items()}
    return stmts, values

def get_eval_drug_ab_combs(data):
    drug_tx = data['prediction'][drug_col]
    antibodies = data['prediction']['Antibody ID']
    drug_ab_combs = zip(drug_tx, antibodies)
    return drug_ab_combs

def run(dec_thresh=0.5, inc_thresh=1.5):
    data = pd.read_data(pd.data_file)
    ab_agents = pd.get_antibody_map(data)

    # If filtering is to be done based on thresholds only,
    # set this to None
    drug_ab_combs = get_eval_drug_ab_combs(data)
    #drug_ab_combs = None

    stmts, values = make_stmts(data, ab_agents, drug_ab_combs=drug_ab_combs,
                               thresh=[dec_thresh, inc_thresh])

    # Now, preassemble the statements to remove duplicates
    pa_dict = preassemble_stmts(stmts)

    with open('data_stmts.pkl', 'wb') as f:
        pickle.dump((pa_dict, values), f, protocol=2)

    return (stmts, values)

if __name__ == '__main__':
    run()
