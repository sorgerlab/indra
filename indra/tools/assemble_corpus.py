from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import sys
try:
    # Python 2
    import cPickle as pickle
except ImportError:
    # Python 3
    import pickle
import logging
from copy import deepcopy, copy
from indra.statements import *
from indra.belief import BeliefEngine
from indra.util import read_unicode_csv
from indra.databases import uniprot_client
from indra.mechlinker import MechLinker
from indra.preassembler import Preassembler
from indra.tools.expand_families import Expander
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler.grounding_mapper import GroundingMapper
from indra.preassembler.grounding_mapper import gm as grounding_map
from indra.preassembler.grounding_mapper import default_agent_map as agent_map
from indra.preassembler.sitemapper import SiteMapper, default_site_map
from copy import deepcopy

logger = logging.getLogger('assemble_corpus')
indra_logger = logging.getLogger('indra').setLevel(logging.DEBUG)

def _filter(kwargs, arg_list):
    return dict(filter(lambda x: x[0] in arg_list, kwargs.items()))

def dump_statements(stmts, fname):
    """Dump a list of statements into a pickle file.

    Parameters
    ----------
    fname : str
        The name of the pickle file to dump statements into.
    """
    if sys.version_info[0] < 3:
        logger.warning('Files pickled in Python 2 may be incompatible with '
                       'Python 3')
    logger.info('Dumping %d statements into %s...' % (len(stmts), fname))
    with open(fname, 'wb') as fh:
        pickle.dump(stmts, fh, protocol=2)

def load_statements(fname, as_dict=False):
    """Load statements from a pickle file.

    Parameters
    ----------
    fname : str
        The name of the pickle file to load statements from.
    as_dict : Optional[bool]
        If True and the pickle file contains a dictionary of statements, it
        is returned as a dictionary. If False, the statements are always
        returned in a list. Default: False

    Returns
    -------
    stmts : list
        A list or dict of statements that were loaded.
    """
    logger.info('Loading %s...' % fname)
    with open(fname, 'rb') as fh:
        # Encoding argument not available in pickle for Python 2
        if sys.version_info[0] < 3:
            stmts = pickle.load(fh)
        # Encoding argument specified here to enable compatibility with
        # pickle files created with Python 2
        else:
            stmts = pickle.load(fh, encoding='latin1')

    if isinstance(stmts, dict):
        if as_dict:
            return stmts
        st = []
        for pmid, st_list in stmts.items():
            st += st_list
        stmts = st
    logger.info('Loaded %d statements' % len(stmts))
    return stmts

def map_grounding(stmts_in, **kwargs):
    """Map grounding using the GroundingMapper.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to map.
    do_rename : Optional[bool]
        If True, Agents are renamed based on their mapped grounding.
    grounding_map : Optional[dict]
        A user supplied grounding map which maps a string to a
        dictionary of database IDs (in the format used by Agents'
        db_refs).
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of mapped statements.
    """
    logger.info('Mapping grounding on %d statements...' % len(stmts_in))
    do_rename = kwargs.get('do_rename')
    gm = kwargs.get('grounding_map', grounding_map)
    if do_rename is None:
        do_rename = True
    gm = GroundingMapper(gm, agent_map)
    stmts_out = gm.map_agents(stmts_in, do_rename=do_rename)
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def map_sequence(stmts_in, **kwargs):
    """Map sequences using the SiteMapper.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to map.
    do_methionine_offset : boolean
        Whether to check for off-by-one errors in site position (possibly)
        attributable to site numbering from mature proteins after
        cleavage of the initial methionine. If True, checks the reference
        sequence for a known modification at 1 site position greater
        than the given one; if there exists such a site, creates the
        mapping. Default is True.
    do_orthology_mapping : boolean
        Whether to check sequence positions for known modification sites
        in mouse or rat sequences (based on PhosphoSitePlus data). If a
        mouse/rat site is found that is linked to a site in the human
        reference sequence, a mapping is created. Default is True.
    do_isoform_mapping : boolean
        Whether to check sequence positions for known modifications
        in other human isoforms of the protein (based on PhosphoSitePlus
        data). If a site is found that is linked to a site in the human
        reference sequence, a mapping is created. Default is True.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of mapped statements.
    """
    logger.info('Mapping sites on %d statements...' % len(stmts_in))
    kwarg_list = ['do_methionine_offset', 'do_orthology_mapping',
                  'do_isoform_mapping']
    sm = SiteMapper(default_site_map)
    valid, mapped = sm.map_sites(stmts_in, **_filter(kwargs, kwarg_list))
    correctly_mapped_stmts = []
    for ms in mapped:
        if all([mm[1] is not None for mm in ms.mapped_mods]):
            correctly_mapped_stmts.append(ms.mapped_stmt)
    stmts_out = valid + correctly_mapped_stmts
    logger.info('%d statements with valid sites' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def run_preassembly(stmts_in, **kwargs):
    """Run preassembly on a list of statements.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to preassemble.
    return_toplevel : Optional[bool]
        If True, only the top-level statements are returned. If False,
        all statements are returned irrespective of level of specificity.
        Default: True
    poolsize : Optional[int]
        The number of worker processes to use to parallelize the
        comparisons performed by the function. If None (default), no
        parallelization is performed. NOTE: Parallelization is only
        available on Python 3.4 and above.
    size_cutoff : Optional[int]
        Groups with size_cutoff or more statements are sent to worker
        processes, while smaller groups are compared in the parent process.
        Default value is 100. Not relevant when parallelization is not
        used.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.
    save_unique : Optional[str]
        The name of a pickle file to save the unique statements into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of preassembled top-level statements.
    """
    dump_pkl_unique = kwargs.get('save_unique')
    be = BeliefEngine()
    pa = Preassembler(hierarchies, stmts_in)
    run_preassembly_duplicate(pa, be, save=dump_pkl_unique)

    dump_pkl = kwargs.get('save')
    return_toplevel = kwargs.get('return_toplevel', True)
    poolsize = kwargs.get('poolsize', None)
    size_cutoff = kwargs.get('size_cutoff', 100)
    options = {'save': dump_pkl, 'return_toplevel': return_toplevel,
               'poolsize': poolsize, 'size_cutoff': size_cutoff}
    stmts_out = run_preassembly_related(pa, be, **options)
    return stmts_out

def run_preassembly_duplicate(preassembler, beliefengine, **kwargs):
    """Run deduplication stage of preassembly on a list of statements.

    Parameters
    ----------
    preassembler : indra.preassembler.Preassembler
        A Preassembler instance
    beliefengine : indra.belief.BeliefEngine
        A BeliefEngine instance
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of unique statements.
    """
    logger.info('Combining duplicates on %d statements...' %
                len(preassembler.stmts))
    dump_pkl = kwargs.get('save')
    stmts_out = preassembler.combine_duplicates()
    beliefengine.set_prior_probs(stmts_out)
    logger.info('%d unique statements' % len(stmts_out))
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def run_preassembly_related(preassembler, beliefengine, **kwargs):
    """Run related stage of preassembly on a list of statements.

    Parameters
    ----------
    preassembler : indra.preassembler.Preassembler
        A Preassembler instance which already has a set of unique statements
        internally.
    beliefengine : indra.belief.BeliefEngine
        A BeliefEngine instance
    return_toplevel : Optional[bool]
        If True, only the top-level statements are returned. If False,
        all statements are returned irrespective of level of specificity.
        Default: True
    poolsize : Optional[int]
        The number of worker processes to use to parallelize the
        comparisons performed by the function. If None (default), no
        parallelization is performed. NOTE: Parallelization is only
        available on Python 3.4 and above.
    size_cutoff : Optional[int]
        Groups with size_cutoff or more statements are sent to worker
        processes, while smaller groups are compared in the parent process.
        Default value is 100. Not relevant when parallelization is not
        used.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of preassembled top-level statements.
    """
    logger.info('Combining related on %d statements...' %
                len(preassembler.unique_stmts))
    return_toplevel = kwargs.get('return_toplevel', True)
    poolsize = kwargs.get('poolsize', None)
    size_cutoff = kwargs.get('size_cutoff', 100)
    stmts_out = preassembler.combine_related(return_toplevel=False,
                                             poolsize=poolsize,
                                             size_cutoff=size_cutoff)
    beliefengine.set_hierarchy_probs(stmts_out)
    stmts_top = filter_top_level(stmts_out)
    if return_toplevel:
        stmts_out = stmts_top
        logger.info('%d top-level statements' % len(stmts_out))
    else:
        logger.info('%d statements out of which %d are top-level' %
                    (len(stmts_out), len(stmts_top)))

    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_by_type(stmts_in, stmt_type, **kwargs):
    """Filter to a given statement type.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    stmt_type : indra.statements.Statement
        The class of the statement type to filter for.
        Example: indra.statements.Modification
    invert : Optional[bool]
        If True, the statements that are not of the given type
        are returned. Default: False
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    invert = kwargs.get('invert', False)
    logger.info('Filtering %d statements for type %s%s...' %
                (len(stmts_in), 'not ' if invert else '',
                 stmt_type.__name__))
    if not invert:
        stmts_out = [st for st in stmts_in if isinstance(st, stmt_type)]
    else:
        stmts_out = [st for st in stmts_in if not isinstance(st, stmt_type)]

    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out


def _agent_is_grounded(agent, score_threshold):
    grounded = True
    db_names = list(set(agent.db_refs.keys()) - set(['TEXT']))
    # If there are no entries at all other than possibly TEXT
    if not db_names:
        grounded = False
    # If there are entries but they point to None / empty values
    if not all([agent.db_refs[db_name] for db_name in db_names]):
        grounded = False
    # If we are looking for scored groundings with a threshold
    if score_threshold:
        for db_name in db_names:
            val = agent.db_refs[db_name]
            # If it's a list with some values, find the
            # highest scoring match and compare to threshold
            if isinstance(val, list) and val:
                high_score = sorted(val, key=lambda x: x[1],
                                    reverse=True)[0][1]
                if high_score < score_threshold:
                    grounded = False
    return grounded


def _remove_bound_conditions(agent, keep_criterion):
    """Removes bound conditions of agent such that keep_criterion is False.

    Parameters
    ----------
    agent: Agent
        The agent whose bound conditions we evaluate
    keep_criterion: function
        Evaluates removal_criterion(a) for each agent a in a bound condition
        and if it evaluates to False, removes a from agent's bound_conditions
    """
    new_bc = []
    for ind in range(len(agent.bound_conditions)):
        if keep_criterion(agent.bound_conditions[ind].agent):
            new_bc.append(agent.bound_conditions[ind])
    agent.bound_conditions = new_bc

def _any_bound_condition_fails_criterion(agent, criterion):
    """Returns True if any bound condition fails to meet the specified
    criterion.

    Parameters
    ----------
    agent: Agent
        The agent whose bound conditions we evaluate
    criterion: function
        Evaluates criterion(a) for each a in a bound condition and returns True
        if any agents fail to meet the criterion.

    Returns
    -------
    any_meets: bool
        True if and only if any of the agents in a bound condition fail to match
        the specified criteria
    """
    bc_agents = [bc.agent for bc in agent.bound_conditions]
    for b in bc_agents:
        if not criterion(b):
            return True
    return False

def filter_grounded_only(stmts_in, **kwargs):
    """Filter to statements that have grounded agents.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    score_threshold : Optional[float]
        If scored groundings are available in a list and the highest score
        if below this threshold, the Statement is filtered out.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.
    remove_bound: Optional[bool]
        If true, removes ungrounded bound conditions from a statement.
        If false (default), filters out statements with ungrounded bound
        conditions.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    remove_bound = kwargs.get('remove_bound', False)

    logger.info('Filtering %d statements for grounded agents...' % 
                len(stmts_in))
    stmts_out = []
    score_threshold = kwargs.get('score_threshold')
    for st in stmts_in:
        grounded = True
        for agent in st.agent_list():
            if agent is not None:
                criterion = lambda x: _agent_is_grounded(x, score_threshold)
                if not criterion(agent):
                    grounded = False
                    break
                if not isinstance(agent, Agent):
                    continue
                if remove_bound:
                    _remove_bound_conditions(agent, criterion)
                elif _any_bound_condition_fails_criterion(agent, criterion):
                    grounded = False
                    break
        if grounded:
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def _agent_is_gene(agent, specific_only):
    """Returns whether an agent is for a gene.

    Parameters
    ----------
    agent: Agent
        The agent to evaluate
    specific_only : Optional[bool]
        If True, only elementary genes/proteins evaluate as genes and families
        will be filtered out. If False, families are also included.

    Returns
    -------
    is_gene: bool
        Whether the agent is a gene
    """
    if not specific_only:
        if not(agent.db_refs.get('HGNC') or \
               agent.db_refs.get('UP') or \
               agent.db_refs.get('FPLX')):
            return False
    else:
        if not(agent.db_refs.get('HGNC') or \
               agent.db_refs.get('UP')):
            return False
    return True

def filter_genes_only(stmts_in, **kwargs):
    """Filter to statements containing genes only.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    specific_only : Optional[bool]
        If True, only elementary genes/proteins will be kept and families
        will be filtered out. If False, families are also included in the
        output. Default: False
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.
    remove_bound: Optional[bool]
        If true, removes bound conditions that are not genes
        If false (default), filters out statements with non-gene bound
        conditions

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """

    if 'remove_bound' in kwargs and kwargs['remove_bound']:
        remove_bound = True
    else:
        remove_bound = False

    specific_only = kwargs.get('specific_only')
    logger.info('Filtering %d statements for ones containing genes only...' % 
                len(stmts_in))
    stmts_out = []
    for st in stmts_in:
        genes_only = True
        for agent in st.agent_list():
            if agent is not None:
                criterion = lambda a: _agent_is_gene(a, specific_only)
                if not criterion(agent):
                    genes_only = False
                    break
                if remove_bound:
                    _remove_bound_conditions(agent, criterion)
                else:
                    if _any_bound_condition_fails_criterion(agent, criterion):
                        genes_only = False
                        break

        if genes_only:
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_belief(stmts_in, belief_cutoff, **kwargs):
    """Filter to statements with belief above a given cutoff.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    belief_cutoff : float
        Only statements with belief above the belief_cutoff will be returned.
        Here 0 < belief_cutoff < 1.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    dump_pkl = kwargs.get('save')
    logger.info('Filtering %d statements to above %f belief' %
                (len(stmts_in), belief_cutoff))
    # The first round of filtering is in the top-level list
    stmts_out = []
    # Now we eliminate supports/supported-by
    for stmt in stmts_in:
        if stmt.belief >= belief_cutoff:
            stmts_out.append(stmt)
        else:
            continue
        supp_by = []
        supp = []
        for st in stmt.supports:
            if st.belief >= belief_cutoff:
                supp.append(st)
        for st in stmt.supported_by:
            if st.belief >= belief_cutoff:
                supp_by.append(st)
        stmt.supports = supp
        stmt.supported_by = supp_by
    logger.info('%d statements after filter...' % len(stmts_out))
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_gene_list(stmts_in, gene_list, policy, allow_families=False,
                     **kwargs):
    """Return statements that contain genes given in a list.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    gene_list : list[str]
        A list of gene symbols to filter for.
    policy : str
        The policy to apply when filtering for the list of genes. "one": keep
        statements that contain at least one of the list of genes and
        possibly others not in the list "all": keep statements that only
        contain genes given in the list
    allow_families : Optional[bool]
        Will include statements involving FamPlex families containing one
        of the genes in the gene list. Default: False
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.
    remove_bound: Optional[str]
        If true, removes bound conditions that are not genes in the list
        If false (default), looks at agents in the bound conditions in addition
        to those participating in the statement directly when applying the
        specified policy.
    invert : Optional[bool]
        If True, the statements that do not match according to the policy
        are returned. Default: False

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    invert = kwargs.get('invert', False)
    remove_bound = kwargs.get('remove_bound', False)

    if policy not in ('one', 'all'):
        logger.error('Policy %s is invalid, not applying filter.' % policy)
    else:
        genes_str = ', '.join(gene_list)
        inv_str = 'not ' if invert else ''
        logger.info(('Filtering %d statements for ones %scontaining "%s" of: '
                     '%s...') % (len(stmts_in), inv_str, policy, genes_str))

    # If we're allowing families, make a list of all FamPlex IDs that
    # contain members of the gene list, and add them to the filter list
    filter_list = copy(gene_list)
    if allow_families:
        for hgnc_name in gene_list:
            gene_uri = hierarchies['entity'].get_uri('HGNC', hgnc_name)
            parents = hierarchies['entity'].get_parents(gene_uri)
            for par_uri in parents:
                ns, id = hierarchies['entity'].ns_id_from_uri(par_uri)
                filter_list.append(id)
    stmts_out = []

    if remove_bound:
        # If requested, remove agents whose names are not in the list from
        # all bound conditions
        if not invert:
            keep_criterion = lambda a: a.name in filter_list
        else:
            keep_criterion = lambda a: a.name not in filter_list

        for st in stmts_in:
            for agent in st.agent_list():
                _remove_bound_conditions(agent, keep_criterion)

    if policy == 'one':
        for st in stmts_in:
            found_gene = False
            if not remove_bound:
                agent_list = st.agent_list_with_bound_condition_agents()
            else:
                agent_list = st.agent_list()
            for agent in agent_list:
                if agent is not None:
                    if agent.name in filter_list:
                        found_gene = True
                        break
            if (found_gene and not invert) or (not found_gene and invert):
                stmts_out.append(st)
    elif policy == 'all':
        for st in stmts_in:
            found_genes = True
            if not remove_bound:
                agent_list = st.agent_list_with_bound_condition_agents()
            else:
                agent_list = st.agent_list()
            for agent in agent_list:
                if agent is not None:
                    if agent.name not in filter_list:
                        found_genes = False
                        break
            if (found_genes and not invert) or (not found_genes and invert):
                stmts_out.append(st)
    else:
        stmts_out = stmts_in

    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out


def filter_human_only(stmts_in, **kwargs):
    """Filter out statements that are grounded, but not to a human gene.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.
    remove_bound: Optional[bool]
        If true, removes all bound conditions that are grounded but not to human
        genes. If false (default), filters out statements with boundary
        conditions that are grounded to non-human genes.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.

    """

    if 'remove_bound' in kwargs and kwargs['remove_bound']:
        remove_bound = True
    else:
        remove_bound = False

    dump_pkl = kwargs.get('save')
    logger.info('Filtering %d statements for human genes only...' %
                len(stmts_in))
    stmts_out = []

    def criterion(agent):
        upid = agent.db_refs.get('UP')
        if upid and not uniprot_client.is_human(upid):
            return False
        else:
            return True


    for st in stmts_in:
        human_genes = True
        for agent in st.agent_list():
            if agent is not None:
                if not criterion(agent):
                    human_genes = False
                    break
                if remove_bound:
                    _remove_bound_conditions(agent, criterion)
                elif _any_bound_condition_fails_criterion(agent, criterion):
                    human_genes = False
                    break
        if human_genes:
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_direct(stmts_in, **kwargs):
    """Filter to statements that are direct interactions

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    def get_is_direct(stmt):
        """Returns true if there is evidence that the statement is a direct
        interaction.

        If any of the evidences associated with the statement
        indicates a direct interatcion then we assume the interaction
        is direct. If there is no evidence for the interaction being indirect
        then we default to direct.
        """
        any_indirect = False
        for ev in stmt.evidence:
            if ev.epistemics.get('direct') is True:
                return True
            elif ev.epistemics.get('direct') is False:
                # This guarantees that we have seen at least
                # some evidence that the statement is indirect
                any_indirect = True
        if any_indirect:
            return False
        return True
    logger.info('Filtering %d statements to direct ones...' % len(stmts_in))
    stmts_out = []
    for st in stmts_in:
        if get_is_direct(st):
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_no_hypothesis(stmts_in, **kwargs):
    """Filter to statements that are not marked as hypothesis in epistemics.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    logger.info('Filtering %d statements to no hypothesis...' % len(stmts_in))
    stmts_out = []
    for st in stmts_in:
        all_hypotheses = True
        ev = None
        for ev in st.evidence:
            if not ev.epistemics.get('hypothesis', False):
                all_hypotheses = False
                break
        if ev is None:
            all_hypotheses = False
        if not all_hypotheses:
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_evidence_source(stmts_in, source_apis, policy='one', **kwargs):
    """Filter to statements that have evidence from a given set of sources.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    source_apis : list[str]
        A list of sources to filter for. Examples: biopax, bel, reach
    policy : Optional[str]
        If 'one', a statement that hase evidence from any of the sources is
        kept. If 'all', only those statements are kept which have evidence
        from all the input sources specified in source_apis.
        If 'none', only those statements are kept that don't have evidence
        from any of the sources specified in source_apis.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    logger.info('Filtering %d statements to evidence source "%s" of: %s...' %
                (len(stmts_in), policy, ', '.join(source_apis)))
    stmts_out = []
    for st in stmts_in:
        sources = set([ev.source_api for ev in st.evidence])
        if policy == 'one':
            if sources.intersection(source_apis):
                stmts_out.append(st)
        if policy == 'all':
            if sources.intersection(source_apis) == set(source_apis):
                stmts_out.append(st)
        if policy == 'none':
            if not sources.intersection(source_apis):
                stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_top_level(stmts_in, **kwargs):
    """Filter to statements that are at the top-level of the hierarchy.

    Here top-level statements correspond to most specific ones.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    logger.info('Filtering %d statements for top-level...' % len(stmts_in))
    stmts_out = [st for st in stmts_in if not st.supports]
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_inconsequential_mods(stmts_in, whitelist=None, **kwargs):
    """Filter out Modifications that modify inconsequential sites

    Inconsequential here means that the site is not mentioned / tested
    in any other statement. In some cases specific sites should be
    preserved, for instance, to be used as readouts in a model.
    In this case, the given sites can be passed in a whitelist.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    whitelist : Optional[dict]
        A whitelist containing agent modification sites whose
        modifications should be preserved even if no other statement
        refers to them. The whitelist parameter is a dictionary in which
        the key is a gene name and the value is a list of tuples of
        (modification_type, residue, position). Example:
        whitelist = {'MAP2K1': [('phosphorylation', 'S', '222')]}
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    if whitelist is None:
        whitelist = {}
    logger.info('Filtering %d statements to remove' % len(stmts_in) +
                ' inconsequential modifications...')
    states_used = whitelist
    for stmt in stmts_in:
        for agent in stmt.agent_list():
            if agent is not None:
                if agent.mods:
                    for mc in agent.mods:
                        mod = (mc.mod_type, mc.residue, mc.position)
                        try:
                            states_used[agent.name].append(mod)
                        except KeyError:
                            states_used[agent.name] = [mod]
    for k, v in states_used.items():
        states_used[k] = list(set(v))
    stmts_out = []
    for stmt in stmts_in:
        skip = False
        if isinstance(stmt, Modification):
            mod_type = modclass_to_modtype[stmt.__class__]
            if isinstance(stmt, RemoveModification):
                mod_type = modtype_to_inverse[mod_type]
            mod = (mod_type, stmt.residue, stmt.position)
            used = states_used.get(stmt.sub.name, [])
            if mod not in used:
                skip = True
        if not skip:
            stmts_out.append(stmt)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_inconsequential_acts(stmts_in, whitelist=None, **kwargs):
    """Filter out Activations that modify inconsequential activities

    Inconsequential here means that the site is not mentioned / tested
    in any other statement. In some cases specific activity types should be
    preserved, for instance, to be used as readouts in a model.
    In this case, the given activities can be passed in a whitelist.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    whitelist : Optional[dict]
        A whitelist containing agent activity types which  should be preserved
        even if no other statement refers to them.
        The whitelist parameter is a dictionary in which
        the key is a gene name and the value is a list of activity types.
        Example: whitelist = {'MAP2K1': ['kinase']}
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    if whitelist is None:
        whitelist = {}
    logger.info('Filtering %d statements to remove' % len(stmts_in) +
                ' inconsequential activations...')
    states_used = whitelist
    for stmt in stmts_in:
        for agent in stmt.agent_list():
            if agent is not None:
                if agent.activity:
                    act = agent.activity.activity_type
                    try:
                        states_used[agent.name].append(act)
                    except KeyError:
                        states_used[agent.name] = [act]
    for k, v in states_used.items():
        states_used[k] = list(set(v))
    stmts_out = []
    for stmt in stmts_in:
        skip = False
        if isinstance(stmt, RegulateActivity):
            used = states_used.get(stmt.obj.name, [])
            if stmt.obj_activity not in used:
                skip = True
        if not skip:
            stmts_out.append(stmt)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def get_unreachable_mods(stmts_in):
    mods_set = {}
    for stmt in stmts_in:
        if isinstance(stmt, Modification):
            mod_type = modclass_to_modtype[stmt.__class__]
            if isinstance(stmt, RemoveModification):
                mod_type = modtype_to_inverse[mod_type]
            mod = (mod_type, stmt.residue, stmt.position)
            if stmt.sub.name not in mods_set:
                mods_set[stmt.sub.name] = set([mod])
            else:
                mods_set[stmt.sub.name].add(mod)
    unreachable_mods = {}
    for stmt in stmts_in:
        for agent in stmt.agent_list():
            if agent is None or not agent.mods:
                continue
            for mc in agent.mods:
                mod = (mc.mod_type, mc.residue, mc.position)
                if mod not in mods_set.get(agent.name, []):
                    msg = '%s not reachable for %s' % (mod, agent.name)
                    logger.warning(msg)
                    if agent.name not in unreachable_mods:
                        unreachable_mods[agent.name] = set([mod])
                    else:
                        unreachable_mods[agent.name].add(mod)

    return unreachable_mods


def filter_mutation_status(stmts_in, mutations, deletions, **kwargs):
    """Filter statements based on existing mutations/deletions

    This filter helps to contextualize a set of statements to a given
    cell type. Given a list of deleted genes, it removes statements that refer
    to these genes. It also takes a list of mutations and removes statements
    that refer to mutations not relevant for the given context.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    mutations : dict
        A dictionary whose keys are gene names, and the values are lists of
        tuples of the form (residue_from, position, residue_to).
        Example: mutations = {'BRAF': [('V', '600', 'E')]}
    deletions : list
        A list of gene names that are deleted.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """

    if 'remove_bound' in kwargs and kwargs['remove_bound']:
        remove_bound = True
    else:
        remove_bound = False

    def criterion(agent):
        if agent is not None and agent.name in deletions:
            return False
        if agent is not None and agent.mutations:
            muts = mutations.get(agent.name, [])
            for mut in agent.mutations:
                mut_tup = (mut.residue_from, mut.position, mut.residue_to)
                if mut_tup not in muts:
                    return False
        return True


    logger.info('Filtering %d statements for mutation status...' %
                len(stmts_in))
    stmts_out = []
    for stmt in stmts_in:
        skip = False
        for agent in stmt.agent_list():
            if not criterion(agent):
                skip = True
                break
            if remove_bound:
                _remove_bound_conditions(agent, criterion)
            elif _any_bound_condition_fails_criterion(agent, criterion):
                skip = True
                break
        if not skip:
            stmts_out.append(stmt)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_enzyme_kinase(stmts_in, **kwargs):
    """Filter Phosphorylations to ones where the enzyme is a known kinase.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    logger.info('Filtering %d statements to remove ' % len(stmts_in) +
                'phosphorylation by non-kinases...')
    path = os.path.dirname(os.path.abspath(__file__))
    kinase_table = read_unicode_csv(path + '/../resources/kinases.tsv',
                                    delimiter='\t')
    gene_names = [lin[1] for lin in list(kinase_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, Phosphorylation):
            if st.enz is not None:
                if st.enz.name in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_mod_nokinase(stmts_in, **kwargs):
    """Filter non-phospho Modifications to ones with a non-kinase enzyme.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    logger.info('Filtering %d statements to remove ' % len(stmts_in) +
                'non-phospho modifications by kinases...')
    path = os.path.dirname(os.path.abspath(__file__))
    kinase_table = read_unicode_csv(path + '/../resources/kinases.tsv',
                                    delimiter='\t')
    gene_names = [lin[1] for lin in list(kinase_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, Modification) and not \
           isinstance(st, Phosphorylation):
            if st.enz is not None:
                if st.enz.name not in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def filter_transcription_factor(stmts_in, **kwargs):
    """Filter out RegulateAmounts where subject is not a transcription factor.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    logger.info('Filtering %d statements to remove ' % len(stmts_in) +
                'amount regulations by non-transcription-factors...')
    path = os.path.dirname(os.path.abspath(__file__))
    tf_table = \
        read_unicode_csv(path + '/../resources/transcription_factors.csv')
    gene_names = [lin[1] for lin in list(tf_table)[1:]]
    stmts_out = []
    for st in stmts_in:
        if isinstance(st, RegulateAmount):
            if st.subj is not None:
                if st.subj.name in gene_names:
                    stmts_out.append(st)
        else:
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out


def filter_uuid_list(stmts_in, uuids, **kwargs):
    """Filter to Statements corresponding to given UUIDs

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to filter.
    uuids : list[str]
        A list of UUIDs to filter for.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of filtered statements.
    """
    logger.info('Filtering %d statements for %d UUID%s...' %
                (len(stmts_in), len(uuids), 's' if len(uuids) > 1 else ''))
    stmts_out = []
    for st in stmts_in:
        if st.uuid in uuids:
            stmts_out.append(st)
    logger.info('%d statements after filter...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out


def expand_families(stmts_in, **kwargs):
    """Expand FamPlex Agents to individual genes.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to expand.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of expanded statements.
    """
    logger.info('Expanding families on %d statements...' % len(stmts_in))
    expander = Expander(hierarchies)
    stmts_out = expander.expand_families(stmts_in)
    logger.info('%d statements after expanding families...' % len(stmts_out))
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def reduce_activities(stmts_in, **kwargs):
    """Reduce the activity types in a list of statements

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to reduce activity types in.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of reduced activity statements.
    """
    logger.info('Reducing activities on %d statements...' % len(stmts_in))
    stmts_out = [deepcopy(st) for st in stmts_in]
    ml = MechLinker(stmts_out)
    ml.gather_explicit_activities()
    ml.reduce_activities()
    stmts_out = ml.statements
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out


def strip_agent_context(stmts_in, **kwargs):
    """Strip any context on agents within each statement.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements whose agent context should be stripped.
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of stripped statements.
    """
    logger.info('Stripping agent context on %d statements...' % len(stmts_in))
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
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out

def dump_stmt_strings(stmts, fname):
    """Save printed statements in a file.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements to save in a text file.
    fname : Optional[str]
        The name of a text file to save the printed statements into.
    """
    with open(fname, 'wb') as fh:
        for st in stmts:
            fh.write(('%s\n' % st).encode('utf-8'))


def rename_db_ref(stmts_in, ns_from, ns_to, **kwargs):
    """Rename an entry in the db_refs of each Agent.

    This is particularly useful when old Statements in pickle files
    need to be updated after a namespace was changed such as
    'BE' to 'FPLX'.

    Parameters
    ----------
    stmts_in : list[indra.statements.Statement]
        A list of statements whose Agents' db_refs need to be changed
    ns_from : str
        The namespace identifier to replace
    ns_to : str
        The namespace identifier to replace to
    save : Optional[str]
        The name of a pickle file to save the results (stmts_out) into.

    Returns
    -------
    stmts_out : list[indra.statements.Statement]
        A list of Statements with Agents' db_refs changed.
    """
    logger.info('Remapping "%s" to "%s" in db_refs on %d statements...' %
                (ns_from, ns_to, len(stmts_in)))
    stmts_out = [deepcopy(st) for st in stmts_in]
    for stmt in stmts_out:
        for agent in stmt.agent_list():
            if agent is not None and ns_from in agent.db_refs:
                agent.db_refs[ns_to] = agent.db_refs.pop(ns_from)
    dump_pkl = kwargs.get('save')
    if dump_pkl:
        dump_statements(stmts_out, dump_pkl)
    return stmts_out


if __name__ == '__main__':
    if len(sys.argv) < 3:
        logger.error('Usage: assemble_corpus.py <pickle_file> <output_folder>')
        sys.exit()
    stmts_fname = sys.argv[1]
    out_folder = sys.argv[2]

    stmts = load_statements(stmts_fname)

    logger.info('All statements: %d' % len(stmts))

    cache_pkl = os.path.join(out_folder, 'mapped_stmts.pkl')
    options = {'save': cache_pkl, 'do_rename': True}
    stmts = map_grounding(stmts, **options)

    cache_pkl = os.path.join(out_folder, 'sequence_valid_stmts.pkl')
    options = {'save': cache_pkl}
    mapped_stmts = map_sequence(stmts, **options)

    be = BeliefEngine()
    pa = Preassembler(hierarchies, mapped_stmts)

    cache_pkl = os.path.join(out_folder, 'unique_stmts.pkl')
    options = {'save': cache_pkl}
    unique_stmts = run_preassembly_duplicate(pa, be, **options)

    cache_pkl = os.path.join(out_folder, 'top_stmts.pkl')
    options = {'save': cache_pkl}
    stmts = run_preassembly_related(pa, be, **options)
