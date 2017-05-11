import process_data as pd
from indra.databases import hgnc_client
from read_phosphosite import read_phosphosite
from indra.statements import Agent, Dephosphorylation, Phosphorylation
from collections import defaultdict
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
import pickle

data_dict = pd.read_data(pd.data_file)
protein_data = data_dict['protein']
abs = pd.get_phos_antibodies(data_dict)
drug_tx = pd.get_single_drug_treatments(data_dict)
drug_targets = pd.get_drug_targets()
ps_data = read_phosphosite('sources/annotated_kinases_v5.csv')

INC_THRESHOLD = 1.2
DEC_THRESHOLD = 0.8
DRUG_COL = 'Sample Description (drug abbre. | dose or time-point)'

stmts = defaultdict(lambda: defaultdict(list))
values = defaultdict(dict)

for tx in drug_tx:
    # Get the drug name
    drug_name = tx.split('|')[0]
    targets = drug_targets[drug_name]
    for target in targets:
        # Create an agent for the drug target
        target_hgnc_id = hgnc_client.get_hgnc_id(target)
        target_up_id = hgnc_client.get_uniprot_id(target_hgnc_id)
        target_agent = Agent(target, db_refs={'HGNC': target_hgnc_id,
                                              'UP': target_up_id})
        for ab in abs:
            try:
                phosphoforms = ps_data[1][ab]
            except KeyError as e:
                print("Error, skipping: %s" % e)
                continue
            drug_tx_data = protein_data[protein_data[DRUG_COL] == tx]
            fold_change = drug_tx_data[ab].values[0]
            if fold_change < DEC_THRESHOLD or fold_change > INC_THRESHOLD:
                values[drug_name][ab] = fold_change
                # Each entry in this list is an Agent with db_refs filled in
                # and associated mod conditions for the phosphorylated sites
                for psf in phosphoforms:
                    psf_agent = Agent(psf.name, db_refs=psf.db_refs)
                    for mod in psf.mods:
                        # Create a Phosphorylation statement corresponding to
                        # this drug/Ab pair
                        if fold_change < DEC_THRESHOLD:
                            stmt = Phosphorylation(target_agent, psf_agent,
                                                   mod.residue, mod.position)
                            stmts[drug_name][ab].append(stmt)
                        else:
                            stmt = Dephosphorylation(target_agent, psf_agent,
                                                     mod.residue, mod.position)
                            stmts[drug_name][ab].append(stmt)

# Now, preassemble the statements to remove duplicates
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

with open('data_stmts.pkl', 'wb') as f:
    pickle.dump((pa_dict, values), f, protocol=2)

