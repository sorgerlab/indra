import process_data as pd
from indra.databases import hgnc_client
from read_phosphosite import read_phosphosite
from indra.statements import Agent, Dephosphorylation, Phosphorylation

data_dict = pd.read_data(pd.data_file)
protein_data = data_dict['protein']
abs = pd.get_phos_antibodies(data_dict)
drug_tx = pd.get_single_drug_treatments(data_dict)
drug_targets = pd.get_drug_targets()
ps_data = read_phosphosite('data/annotated_kinases_v4.csv')

INC_THRESHOLD = 1.5
DEC_THRESHOLD = 0.5
DRUG_COL = 'Sample Description (drug abbre. | dose or time-point)'

stmts = []


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
            except KeyError:
                continue
            # Each entry in this list is an Agent with db_refs filled in
            # and associated mod conditions for the phosphorylated sites
            for psf in phosphoforms:
                psf_agent = Agent(psf.name, db_refs=psf.db_refs)
                for mod in psf.mods:
                    # Create a Phosphorylation statement corresponding to this
                    # drug/Ab pair
                    drug_tx_data = protein_data[protein_data[DRUG_COL] == tx]
                    fold_change = drug_tx_data[ab].values[0]
                    if fold_change < DEC_THRESHOLD:
                        stmt = Phosphorylation(target_agent, psf_agent,
                                               mod.residue, mod.position)
                        stmts.append(stmt)
                    elif fold_change > INC_THRESHOLD:
                        stmt = Dephosphorylation(target_agent, psf_agent,
                                                 mod.residue, mod.position)
                        stmts.append(stmt)

            # Get the target of the drug
            # Determine if the ab measurement went up or down
            # Get the modified protein states associated with the ab

    # Get the data for the drug row
    # Iterate
