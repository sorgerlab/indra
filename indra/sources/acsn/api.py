__all__ = ['process_from_web', 'process_df']

import os
import pandas as pd
from collections import defaultdict
from .processor import AcsnProcessor
from indra.ontology.bio import bio_ontology

ACSN_URL = 'https://acsn.curie.fr/ACSN2/downloads/'
ACSN_RELATIONS_URL = ACSN_URL + \
                     'ACSN2_binary_relations_between_proteins_with_PMID.txt'
ACSN_CORRESPONDENCE = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                                   'ACSN2_HUGO_Correspondence.gmt')


def transform_gmt(gmt):
    # Convert the GMT file into a dictionary
    acsn_gmt_dict = defaultdict(set)
    with open(gmt, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            line = line.rstrip('\t')
            interactors = line.split('\t')
            for ag2 in interactors[2:]:
                acsn_gmt_dict[(interactors[0])].add(ag2)
    return acsn_gmt_dict


def process_from_web():
    relations_df = pd.read_csv(ACSN_RELATIONS_URL, sep='\t')
    correspondence_dict = transform_gmt(ACSN_CORRESPONDENCE)
    return process_df(relations_df, correspondence_dict)


def process_df(relations_df, correspondence_dict):
    fplx_lookup = get_famplex_lookup()
    ap = AcsnProcessor(relations_df, correspondence_dict,
                       fplx_lookup)
    ap.extract_statements()
    return ap


def get_famplex_lookup():
    fplx_lookup = {}
    bio_ontology.initialize()
    for node in bio_ontology.nodes:
        ns, id = bio_ontology.get_ns_id(node)
        if ns == 'FPLX':
            children = bio_ontology.get_children(ns, id)
            hgnc_children = [bio_ontology.get_name(*c)
                             for c in children if c[0] == 'HGNC']
            fplx_lookup[tuple(sorted(hgnc_children))] = id
    return fplx_lookup
