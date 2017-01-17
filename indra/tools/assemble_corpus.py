from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import sys
import time
try:
    # Python 2
    import cPickle as pickle
except ImportError:
    # Python 3
    import pickle
import logging
from copy import deepcopy
from indra.statements import *
from indra.belief import BeliefEngine
from indra.databases import uniprot_client
from indra.preassembler import Preassembler
from indra.tools.expand_families import Expander
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler.grounding_mapper import GroundingMapper
from indra.preassembler.grounding_mapper import gm as grounding_map
from indra.preassembler.sitemapper import SiteMapper, default_site_map

logger = logging.getLogger('assemble_corpus')
indra_logger = logging.getLogger('indra').setLevel(logging.DEBUG)

def dump_statements(stmts, fname):
    """Dump a list of statements into a pickle file."""
    logger.info('Dumping into %s' % fname)
    with open(fname, 'wb') as fh:
        pickle.dump(stmts, fh, protocol=2)

def load_statements(fname):
    """Load a list of statements from a pickle file."""
    logger.info('Loading %s' % fname)
    with open(fname, 'rb') as fh:
        stmts = pickle.load(fh)
    if isinstance(stmts, dict):
        st = []
        for pmid, st_list in stmts.items():
            st += st_list
        stmts = st
    return stmts

def map_grounding(stmts_in, **kwargs):
    """Map grounding using the GroundingMapper."""
    logger.info('Mapping grounding...')
    load_pkl = kwargs.get('load_pkl')
    dump_pkl = kwargs.get('dump_pkl')
    do_rename = kwargs.get('do_rename')
    if do_rename is None:
        do_rename = True
    if load_pkl:
        stmts_out = load_statements(load_pkl)
        return stmts_out
    gm = GroundingMapper(grounding_map)
    stmts_out = gm.map_agents(stmts_in, do_rename=do_rename)
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def map_sequence(stmts_in, **kwargs):
    """Map sequences using the SiteMapper."""
    logger.info('Mapping sites...')
    load_pkl = kwargs.get('load_pkl')
    dump_pkl = kwargs.get('dump_pkl')
    if load_pkl:
        stmts_out = load_statements(load_pkl)
    else:
        sm = SiteMapper(default_site_map)
        stmts_out, _ = sm.map_sites(stmts_in)
        if dump_pkl:
            dump_statements(stmts_out, dump_pkl)
    logger.info('Statements with valid sites: %d' % len(stmts_out))
    return stmts_out

def run_preassembly(stmts_in, **kwargs):
    """Run preassembly."""
    dump_pkl = kwargs.get('dump_pkl')
    be = BeliefEngine()
    pa = Preassembler(hierarchies, stmts_in)

    options = {'preassembler': pa, 'beliefengine': be}
    unique_stmts = run_preassembly_duplicate(stmts_in, **options)

    options = {'dump_pkl': dump_pkl, 'preassembler': pa, 'beliefengine': be}
    stmts_out = run_preassembly_related(unique_stmts, **options)

    return stmts_out

def run_preassembly_duplicate(stmts_in, **kwargs):
    """Run deduplication stage of preassembly."""
    logger.info('Combining duplicates...')
    load_pkl = kwargs.get('load_pkl')
    dump_pkl = kwargs.get('dump_pkl')
    pa = kwargs.get('preassembler')
    be = kwargs.get('beliefengine')
    if load_pkl:
        stmts_out = load_statements(load_pkl)
    else:
        stmts_out = pa.combine_duplicates()
        if dump_pkl:
            dump_statements(stmts_out, dump_pkl)
    logger.info('Unique statements: %d' % len(stmts_out))
    be.set_prior_probs(stmts_out)
    return stmts_out

def run_preassembly_related(stmts_in, **kwargs):
    """Run related stage of preassembly."""
    logger.info('Combining related...')
    load_pkl = kwargs.get('load_pkl')
    dump_pkl = kwargs.get('dump_pkl')
    pa = kwargs.get('preassembler')
    be = kwargs.get('beliefengine')
    if load_pkl:
        stmts_out = load_statements(load_pkl)
    else:
        stmts_out = pa.combine_related()
        if dump_pkl:
            dump_statements(stmts_out, dump_pkl)
    logger.info('Top-level statements: %d' % len(stmts_out))
    be.set_hierarchy_probs(stmts_out)
    return stmts_out

def filter_by_type(stmts_in, stmt_type):
    """Filter to a given statement type."""
    stmts_out = [st for st in stmts_in if isinstance(st, stmt_type)]
    return stmts_out

def filter_grounded_only(stmts_in, **kwargs):
    """Filter to statements that have grounded agents."""
    load_pkl = kwargs.get('load_pkl')
    dump_pkl = kwargs.get('dump_pkl')
    logger.info('Filtering %d statements for grounded agents...' % 
                len(stmts_in))
    if load_pkl:
        stmts_out = load_statements(load_pkl)
        return stmts_out
    stmts_out = []
    for st in stmts_in:
        grounded = True
        for agent in st.agent_list():
            if agent is not None:
                if (len(agent.db_refs) == 1) and agent.db_refs.get('TEXT'):
                    grounded = False
                    break
        if grounded:
            stmts_out.append(st)
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_genes_only(stmts_in, **kwargs):
    """Filter to statements containing genes only."""
    load_pkl = kwargs.get('load_pkl')
    dump_pkl = kwargs.get('dump_pkl')
    specific_only = kwargs.get('specific_only')
    logger.info('Filtering %d statements for ones containing genes only...' % 
                len(stmts_in))
    if load_pkl:
        stmts_out = load_statements(load_pkl)
        return stmts_out
    stmts_out = []
    for st in stmts_in:
        genes_only = True
        for agent in st.agent_list():
            if agent is not None:
                if not specific_only:
                    if not(agent.db_refs.get('HGNC') or \
                        agent.db_refs.get('UP') or \
                        agent.db_refs.get('BE')):
                        genes_only = False
                        break
                else:
                    if not(agent.db_refs.get('HGNC') or \
                        agent.db_refs.get('UP')):
                        genes_only = False
                        break
        if genes_only:
            stmts_out.append(st)
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_belief(stmts_in, belief_cutoff, **kwargs):
    """Filter to statements with belief about a cutoff."""
    stmts_out = [s for s in stmts_in if s.belief >= belief_cutoff]
    return stmts_out

def expand_families(stmts_in, **kwargs):
    """Expand Bioentities Agents to individual genes."""
    expander = Expander(hierarchies)
    stmts_out = expander.expand_families(stmts_in)
    return stmts_out

def filter_gene_list(stmts_in, gene_list, policy, **kwargs):
    """Return statements that contain genes given in a list.

    gene_list : list[str]
        A list of gene symbols to filter for.
    policy : str
        The policy to apply when filtering for the list of genes.
        'one': keep statements that contain at least one of the
               list of genes and possibly others not in the list
        'all': keep statements that only contain genes given in the list
    """
    load_pkl = kwargs.get('load_pkl')
    dump_pkl = kwargs.get('dump_pkl')
    if policy not in ('one', 'all'):
        logger.error('Policy %s is invalid, not applying filter.' % policy)
    genes_str = ', '.join(gene_list)
    logger.info('Filtering %d statements for ones containing: %s...' %
                (len(stmts_in), genes_str))
    if load_pkl:
        stmts_out = load_statements(load_pkl)
        return stmts_out
    stmts_out = []
    if policy == 'one':
        for st in stmts_in:
            found_gene = False
            for agent in st.agent_list():
                if agent is not None:
                    if agent.name in gene_list:
                        found_gene = True
                        break
            if found_gene:
                stmts_out.append(st)
    elif policy == 'all':
        for st in stmts_in:
            found_genes = True
            for agent in st.agent_list():
                if agent is not None:
                    if agent.name not in gene_list:
                        found_genes = False
                        break
            if found_genes:
                stmts_out.append(st)
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out


def filter_human_only(stmts_in, **kwargs):
    """Filter out statements that are not grounded to human genes."""
    load_pkl = kwargs.get('load_pkl')
    dump_pkl = kwargs.get('dump_pkl')
    logger.info('Filtering %d statements for human genes only...' % 
                len(stmts_in))
    if load_pkl:
        stmts_out = load_statements(load_pkl)
        return stmts_out
    stmts_out = []
    for st in stmts_in:
        human_genes = True
        for agent in st.agent_list():
            if agent is not None:
                upid = agent.db_refs.get('UP')
                if upid and not uniprot_client.is_human(upid):
                    human_genes = False
                    break
        if human_genes:
            stmts_out.append(st)
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def strip_agent_context(stmts_in, **kwargs):
    """Strip any context on agents within each statement."""
    stmts_out = []
    for st in stmts_in:
        new_st = deepcopy(st)
        for agent in new_st.agent_list():
            if agent is None:
                continue
            agent.mods = []
            agent.mutations = []
            agent.activity = None
            agent.location = None
            agent.bound_conditions = []
        stmts_out.append(new_st)
    return stmts_out

def dump_stmt_strings(stmts):
    """Save printed statements in a file."""
    with open('%s.txt' % len(stmts), 'wb') as fh:
        for st in stmts:
            fh.write(('%s\n' % st).encode('utf-8'))

if __name__ == '__main__':
    if len(sys.argv) < 3:
        logger.error('Usage: assemble_corpus.py <pickle_file> <output_folder>')
        sys.exit()
    stmts_fname = sys.argv[1]
    out_folder = sys.argv[2]

    stmts = load_statements(stmts_fname)

    logger.info('All statements: %d' % len(stmts))

    cache_pkl = os.path.join(out_folder, 'mapped_stmts.pkl')
    options = {'dump_pkl': cache_pkl,
               'do_rename': True}
    stmts = map_grounding(stmts, **options)

    cache_pkl = os.path.join(out_folder, 'sequence_valid_stmts.pkl')
    options = {'dump_pkl': cache_pkl}
    mapped_stmts = map_sequence(stmts, **options)

    be = BeliefEngine()
    pa = Preassembler(hierarchies, mapped_stmts)

    cache_pkl = os.path.join(out_folder, 'unique_stmts.pkl')
    options = {'dump_pkl': cache_pkl, 'preassembler': pa, 'beliefengine': be}
    unique_stmts = run_preassembly_duplicate(mapped_stmts, **options)

    cache_pkl = os.path.join(out_folder, 'top_stmts.pkl')
    options = {'dump_pkl': cache_pkl, 'preassembler': pa, 'beliefengine': be}
    stmts = run_preassembly_related(unique_stmts, **options)
