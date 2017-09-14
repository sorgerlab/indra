import os
import io
import sys
import json
import time
import pickle
import numpy as np
from random import shuffle
from matplotlib import pyplot as plt
from indra.preassembler import grounding_mapper as gm

import ndex
from indra.sources import ndex_cx
from indra.databases import hgnc_client
import indra.tools.assemble_corpus as ac
from indra.assemblers import CxAssembler
from indra.literature.pubmed_client import get_ids_for_gene
from indra.util import _require_python3
from indra.tools.gene_network import GeneNetwork
from indra.statements import IncreaseAmount

def build_prior(genes, file_prefix):
    gn = GeneNetwork(genes, file_prefix)
    stmts = gn.get_statements(filter=False)
    ac.dump_statements(stmts, '%s.pkl' % file_prefix)
    return stmts


def get_pmids(gene_names):
    pmids = []
    for gene_name in gene_names:
        pm = get_ids_for_gene(gene_name)
        pmids += pm
        print('%s: %d PMIDs' % (gene_name, len(pm)))
    return pmids


def save_pmids_for_reading(pmids, fname):
    shuffle(pmids)
    with open(fname, 'wt') as fh:
        for pmid in pmids:
            fh.write('%s\n' % pmid)


def run_assembly(stmts, filename):
    stmts = ac.map_grounding(stmts)
    stmts = ac.filter_grounded_only(stmts)
    stmts = ac.filter_human_only(stmts)
    #stmts = ac.expand_families(stmts)
    stmts = ac.filter_gene_list(stmts, gene_names, 'one', allow_families=True)
    #stmts = ac.filter_gene_list(stmts, gene_names, 'all', allow_families=True)
    stmts = ac.map_sequence(stmts)
    stmts = ac.run_preassembly(stmts, return_toplevel=False, poolsize=4)
    ac.dump_statements(stmts, filename)
    return stmts


def filter(stmts, cutoff, filename):
    stmts = ac.filter_belief(stmts, cutoff)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_direct(stmts)
    #stmts = ac.filter_enzyme_kinase(stmts)
    ac.dump_statements(stmts, filename)
    return stmts


def assemble_cx(stmts, save_file):
    cxa = CxAssembler(stmts)
    cxa.make_model(add_indra_json=False)
    cxa.save_model(save_file)
    return cxa


def plot_belief_scores(stmts):
    scores = np.array([s.belief for s in stmts])
    plt.hist(scores)


def upload_to_ndex(cx_str, ndex_cred, network_id):
    server = 'http://public.ndexbio.org'
    username = ndex_cred.get('username')
    password = ndex_cred.get('password')
    nd = ndex.client.Ndex(server, username, password)
    #network_id = ndex_cred.get('network')

    try:
        print('Getting network summary...')
        summary = nd.get_network_summary(network_id)
    except Exception as e:
        print('Could not get NDEx network summary.')
        print(e)
        return

    # Update network content
    try:
        print('Updating network...')
        cx_stream = io.BytesIO(cx_str.encode('utf-8'))
        nd.update_cx_network(cx_stream, network_id)
    except Exception as e:
        print('Could not update NDEx network.')
        print(e)
        return

    # Update network profile
    ver_str = summary.get('version')
    new_ver = _increment_ndex_ver(ver_str)
    profile = {'name': summary.get('name'),
               'description': summary.get('description'),
               'version': new_ver,
               }
    print('Updating NDEx network (%s) profile to %s' %
                (network_id, profile))
    profile_retries = 5
    for i in range(profile_retries):
        try:
            time.sleep(5)
            nd.update_network_profile(network_id, profile)
            break
        except Exception as e:
            print('Could not update NDEx network profile.')
            print(e)

    # Update network style
    import ndex.beta.toolbox as toolbox
    template_uuid = "ea4ea3b7-6903-11e7-961c-0ac135e8bacf"

    d_edge_types = ["Activation", "Inhibition",
                    "Modification", "SelfModification",
                    "Gap", "Gef", "IncreaseAmount",
                    "DecreaseAmount"]

    source_network = ndex.networkn.NdexGraph(server=server, username=username,
                                             password=password,
                                             uuid=network_id)

    toolbox.apply_template(source_network, template_uuid, server=server,
                           username=username, password=password)

    source_network.update_to(network_id, server=server, username=username,
                             password=password)

def _increment_ndex_ver(ver_str):
    if not ver_str:
        new_ver = '1.0'
    else:
        major_ver, minor_ver = ver_str.split('.')
        new_minor_ver = str(int(minor_ver) + 1)
        new_ver = major_ver + '.' + new_minor_ver
    return new_ver

if __name__ == '__main__':
    # Load NDEx credentials
    with open('ndex_cred.json', 'rt') as f:
        ndex_cred = json.load(f)

    # Get the network
    ncp = ndex_cx.process_ndex_network('55cfccd4-966e-11e7-a10d-0ac135e8bacf',
                                       username=ndex_cred['username'],
                                       password=ndex_cred['password'],
                                       require_grounding=False)

    gnd_map_ext = {'PALB2_wt_eto': {'HGNC': 'PALB2'},
                   'PALB2_wt': {'HGNC': 'PALB2'},
                   'CSDA': {'HGNC': 'YBX3'},
                   'COBRA1': {'HGNC': 'NELFB'},
                   'SRPR': {'HGNC': 'SRPRA'},
                   'TOMM70A': {'HGNC': 'TOMM70'}
                   }

    gm.default_grounding_map.update(gnd_map_ext)
    gmapper = gm.GroundingMapper(gm.default_grounding_map)
    ncp_stmts = gmapper.map_agents(ncp.statements)

    gene_names = [hgnc_client.get_hgnc_name(ag.db_refs['HGNC'])
                  for stmt in ncp_stmts for ag in stmt.agent_list()]

    """
    # Get PMIDs for reading
    entrez_pmids = get_pmids(gene_names)
    network_pmids = ncp.get_pmids()
    pmids = list(set(entrez_pmids + network_pmids))
    save_pmids_for_reading(pmids, 'pmids.txt')
    """

    # Build the model
    prior_stmts = build_prior(gene_names, 'palb2_prior')
    #with open('palb2_prior.pkl', 'rb') as f:
    #    prior_stmts = pickle.load(f)
    #    print("Loaded %d prior stmts" % len(prior_stmts))
    reach_stmts = ac.load_statements('reach_stmts.pkl')
    stmts = ncp_stmts + reach_stmts + prior_stmts
    stmts = run_assembly(stmts, 'unfiltered_assembled_stmts.pkl')

    # Filter the statements at different levels
    ids_cutoffs = (('c793e3d1-97f1-11e7-a10d-0ac135e8bacf', 0.90),
                   ('c816a864-97f1-11e7-a10d-0ac135e8bacf', 0.95),
                   ('c83b6e77-97f1-11e7-a10d-0ac135e8bacf', 0.99))

    for net_id, cutoff in ids_cutoffs:
        stmts_filt = filter(stmts, cutoff, 'palb2_stmts_%.2f.pkl' % cutoff)
        cxa = assemble_cx(stmts_filt, 'palb2_%.2f.cx' % cutoff)
        cx_str = cxa.print_cx()
        upload_to_ndex(cx_str, ndex_cred, net_id)
