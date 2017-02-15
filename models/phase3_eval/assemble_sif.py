from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import itertools
from copy import deepcopy
from os.path import join as pjoin
from indra.assemblers import SifAssembler, CxAssembler
import indra.tools.assemble_corpus as ac
from indra.statements import *

import process_data

def assemble_sif(stmts, data, out_file):
    """Return an assembled SIF."""
    # Filter for high-belief statements
    stmts = ac.filter_belief(stmts, 0.99)
    stmts = ac.filter_top_level(stmts)
    # Filter for Activation / Inhibition
    stmts_act = ac.filter_by_type(stmts, Activation)
    stmts_inact = ac.filter_by_type(stmts, Inhibition)
    stmts = stmts_act + stmts_inact
    # Get Ras227 and filter statments
    ras_genes = process_data.get_ras227_genes()
    #ras_genes = [x for x in ras_genes if x not in ['YAP1']]
    stmts = ac.filter_gene_list(stmts, ras_genes, 'all')
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
    stmts = filter_ab_edges(stmts, 'all')
    # Get a list of the AB names that end up being covered in the prior network
    # This is important because other ABs will need to be taken out of the
    # MIDAS file to work.
    def get_ab_names(st):
        prior_abs = set()
        for s in st:
            for a in s.agent_list():
                if a is not None:
                    if a.name.find('_p') != -1:
                        prior_abs.add(a.name)
        return sorted(list(prior_abs))
    pkn_abs = get_ab_names(stmts)
    def get_drug_names(st):
        prior_drugs = set()
        for s in st:
            for a in s.agent_list():
                if a is not None:
                    if a.name.find('Drugs') != -1:
                        prior_drugs.add(a.name.split(':')[0])
        return sorted(list(prior_drugs))
    pkn_drugs = get_drug_names(stmts)
    print('Boolean PKN contains these antibodies: %s' % ', '.join(pkn_abs))
    # Because of a bug in CNO,
    # node names containing AND need to be replaced
    # node names containing - need to be replaced
    # node names starting in a digit need to be replaced
    # must happen before SIF assembly, but not sooner as that will drop
    # names from the MIDAS file
    def rename_nodes(st):
        for s in st:
            for a in s.agent_list():
                if a is not None:
                    if a.name.find('AND') != -1:
                        a.name = a.name.replace('AND', 'A_ND')
                    if a.name.find('-') != -1:
                        a.name = a.name.replace('-', '_')
                    if a.name[0].isdigit():
                        a.name = 'abc_' + a.name
    rename_nodes(stmts)
    # Make the SIF model
    sa = SifAssembler(stmts)
    sa.make_model(use_name_as_key=True)
    sif_str = sa.print_model()
    # assemble and dump a cx of the sif
    ca = CxAssembler()
    ca.add_statements(stmts)
    model = ca.make_model()
    ca.save_model('sif.cx')
    with open(out_file, 'wb') as fh:
        fh.write(sif_str.encode('utf-8'))
    # Make the MIDAS data file used for training the model
    midas_data = process_data.get_midas_data(data, pkn_abs, pkn_drugs)
    return sif_str

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
