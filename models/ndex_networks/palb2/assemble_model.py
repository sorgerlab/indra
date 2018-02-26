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

from indra.sources import ndex_cx
from indra.databases import hgnc_client, ndex_client
import indra.tools.assemble_corpus as ac
from indra.assemblers import CxAssembler
from indra.literature.pubmed_client import get_ids_for_gene
from indra.util import _require_python3
from indra.tools.gene_network import GeneNetwork
from indra.statements import IncreaseAmount

def build_prior(genes, file_prefix):
    gn = GeneNetwork(genes, file_prefix)
    #stmts = gn.get_statements(filter=False)
    stmts = gn.get_biopax_stmts(filter=False)
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


if __name__ == '__main__':
    # Load NDEx credentials
    with open('ndex_cred.json', 'rt') as f:
        ndex_cred = json.load(f)

    # Get the network
    ncp = ndex_cx.process_ndex_network('55cfccd4-966e-11e7-a10d-0ac135e8bacf',
                                       username=ndex_cred['user'],
                                       password=ndex_cred['password'],
                                       require_grounding=False)
    # Add grounding entries for ungrounded nodes/noncanonical gene names
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
        ndex_client.update_network(cx_str, net_id, ndex_cred)
