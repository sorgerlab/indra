__all__ = ['get_statements', 'get_statements_for_paper',
           'get_statements_by_hash', 'submit_curation']

from indra.util import clockit
from indra.statements import stmts_from_json, Complex, SelfModification, \
    ActiveForm, Translocation, Conversion

from indra.sources.indra_db_rest.processor import IndraDBRestSearchProcessor, \
    IndraDBRestHashProcessor
from indra.sources.indra_db_rest.util import submit_statement_request, \
    make_db_rest_request, get_url_base


@clockit
def get_statements(subject=None, object=None, agents=None, stmt_type=None,
                   use_exact_type=False, persist=True, timeout=None,
                   simple_response=False, ev_limit=10, best_first=True, tries=2,
                   max_stmts=None):
    """Get a processor for the INDRA DB web API matching given agents and type.

    There are two types of responses available. You can just get a list of
    INDRA Statements, or you can get an IndraDBRestProcessor object, which allow
    Statements to be loaded in a background thread, providing a sample of the
    best* content available promptly in the sample_statements attribute, and
    populates the statements attribute when the paged load is complete.

    The latter should be used in all new code, and where convenient the prior
    should be converted to use the processor, as this option may be removed in
    the future.

    * In the sense of having the most supporting evidence.

    Parameters
    ----------
    subject/object : str
        Optionally specify the subject and/or object of the statements in
        you wish to get from the database. By default, the namespace is assumed
        to be HGNC gene names, however you may specify another namespace by
        including `@<namespace>` at the end of the name string. For example, if
        you want to specify an agent by chebi, you could use `CHEBI:6801@CHEBI`,
        or if you wanted to use the HGNC id, you could use `6871@HGNC`.
    agents : list[str]
        A list of agents, specified in the same manner as subject and object,
        but without specifying their grammatical position.
    stmt_type : str
        Specify the types of interactions you are interested in, as indicated
        by the sub-classes of INDRA's Statements. This argument is *not* case
        sensitive. If the statement class given has sub-classes
        (e.g. RegulateAmount has IncreaseAmount and DecreaseAmount), then both
        the class itself, and its subclasses, will be queried, by default. If
        you do not want this behavior, set use_exact_type=True. Note that if
        max_stmts is set, it is possible only the exact statement type will
        be returned, as this is the first searched. The processor then cycles
        through the types, getting a page of results for each type and adding it
        to the quota, until the max number of statements is reached.
    use_exact_type : bool
        If stmt_type is given, and you only want to search for that specific
        statement type, set this to True. Default is False.
    persist : bool
        Default is True. When False, if a query comes back limited (not all
        results returned), just give up and pass along what was returned.
        Otherwise, make further queries to get the rest of the data (which may
        take some time).
    timeout : positive int or None
        If an int, block until the work is done and statements are retrieved, or
        until the timeout has expired, in which case the results so far will be
        returned in the response object, and further results will be added in
        a separate thread as they become available. If simple_response is True,
        all statements available will be returned. Otherwise (if None), block
        indefinitely until all statements are retrieved. Default is None.
    simple_response : bool
        If True, a simple list of statements is returned (thus block should also
        be True). If block is False, only the original sample will be returned
        (as though persist was False), until the statements are done loading, in
        which case the rest should appear in the list. This behavior is not
        encouraged. Default is False (which breaks backwards compatibility with
        usage of INDRA versions from before 1/22/2019). WE ENCOURAGE ALL NEW
        USE-CASES TO USE THE PROCESSOR, AS THIS FEATURE MAY BE REMOVED AT A
        LATER DATE.
    ev_limit : int or None
        Limit the amount of evidence returned per Statement. Default is 10.
    best_first : bool
        If True, the preassembled statements will be sorted by the amount of
        evidence they have, and those with the most evidence will be
        prioritized. When using `max_stmts`, this means you will get the "best"
        statements. If False, statements will be queried in arbitrary order.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can also
        help gracefully handle an unreliable connection, if you're willing to
        wait. Default is 2.
    max_stmts : int or None
        Select the maximum number of statements to return. When set less than
        1000 the effect is much the same as setting persist to false, and will
        guarantee a faster response. Default is None.

    Returns
    -------
    processor : :py:class:`IndraDBRestSearchProcessor`
        An instance of the IndraDBRestProcessor, which has an attribute
        `statements` which will be populated when the query/queries are done.
        This is the default behavior, and is encouraged in all future cases,
        however a simple list of statements may be returned using the
        `simple_response` option described above.
    """
    processor = IndraDBRestSearchProcessor(subject, object, agents, stmt_type,
                                           use_exact_type, persist, timeout,
                                           ev_limit, best_first, tries, max_stmts)

    # Format the result appropriately.
    if simple_response:
        ret = processor.statements
    else:
        ret = processor
    return ret


@clockit
def get_statements_by_hash(hash_list, ev_limit=100, best_first=True, tries=2,
                           simple_response=False):
    """Get fully formed statements from a list of hashes.

    Parameters
    ----------
    hash_list : list[int or str]
        A list of statement hashes.
    ev_limit : int or None
        Limit the amount of evidence returned per Statement. Default is 100.
    best_first : bool
        If True, the preassembled statements will be sorted by the amount of
        evidence they have, and those with the most evidence will be
        prioritized. When using `max_stmts`, this means you will get the "best"
        statements. If False, statements will be queried in arbitrary order.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can
        also help gracefully handle an unreliable connection, if you're
        willing to wait. Default is 2.
    simple_response : bool
        If True, a simple list of statements is returned (thus block should also
        be True). If block is False, only the original sample will be returned
        (as though persist was False), until the statements are done loading, in
        which case the rest should appear in the list. This behavior is not
        encouraged. Default is False (which breaks backwards compatibility with
        usage of INDRA versions from before 9/19/2019). WE ENCOURAGE ALL NEW
        USE-CASES TO USE THE PROCESSOR, AS THIS FEATURE MAY BE REMOVED AT A
        LATER DATE.
    """
    processor = IndraDBRestHashProcessor(hash_list, ev_limit=ev_limit,
                                         best_first=best_first, tries=tries)

    if simple_response:
        ret = processor.statements
    else:
        ret = processor
    return ret


@clockit
def get_statements_for_paper(ids, ev_limit=10, best_first=True, tries=2,
                             max_stmts=None):
    """Get the set of raw Statements extracted from a paper given by the id.

    Parameters
    ----------
    ids : list[(<id type>, <id value>)]
        A list of tuples with ids and their type. The type can be any one of
        'pmid', 'pmcid', 'doi', 'pii', 'manuscript id', or 'trid', which is the
        primary key id of the text references in the database.
    ev_limit : int or None
        Limit the amount of evidence returned per Statement. Default is 10.
    best_first : bool
        If True, the preassembled statements will be sorted by the amount of
        evidence they have, and those with the most evidence will be
        prioritized. When using `max_stmts`, this means you will get the "best"
        statements. If False, statements will be queried in arbitrary order.
    tries : int > 0
        Set the number of times to try the query. The database often caches
        results, so if a query times out the first time, trying again after a
        timeout will often succeed fast enough to avoid a timeout. This can also
        help gracefully handle an unreliable connection, if you're willing to
        wait. Default is 2.
    max_stmts : int or None
        Select a maximum number of statements to be returned. Default is None.

    Returns
    -------
    stmts : list[:py:class:`indra.statements.Statement`]
        A list of INDRA Statement instances.
    """
    id_l = [{'id': id_val, 'type': id_type} for id_type, id_val in ids]
    resp = submit_statement_request('post', 'from_papers', data={'ids': id_l},
                                    ev_limit=ev_limit, best_first=best_first,
                                    tries=tries, max_stmts=max_stmts)
    stmts_json = resp.json()['statements']
    return stmts_from_json(stmts_json.values())


def submit_curation(hash_val, tag, curator, text=None,
                    source='indra_rest_client', ev_hash=None, is_test=False):
    """Submit a curation for the given statement at the relevant level.

    Parameters
    ----------
    hash_val : int
        The hash corresponding to the statement.
    tag : str
        A very short phrase categorizing the error or type of curation,
        e.g. "grounding" for a grounding error, or "correct" if you are
        marking a statement as correct.
    curator : str
        The name or identifier for the curator.
    text : str
        A brief description of the problem.
    source : str
        The name of the access point through which the curation was performed.
        The default is 'direct_client', meaning this function was used
        directly. Any higher-level application should identify itself here.
    ev_hash : int
        A hash of the sentence and other evidence information. Elsewhere
        referred to as `source_hash`.
    is_test : bool
        Used in testing. If True, no curation will actually be added to the
        database.
    """
    data = {'tag': tag, 'text': text, 'curator': curator, 'source': source,
            'ev_hash': ev_hash}
    url = 'curation/submit/%s' % hash_val
    if is_test:
        qstr = '?test'
    else:
        qstr = ''
    return make_db_rest_request('post', url, qstr, data=data)


def get_statement_queries(stmts, **params):
    """Get queries used to search based on a statement.

    In addition to the stmts, you can enter any parameters standard to the
    query. See https://github.com/indralab/indra_db/rest_api for a full list.

    Parameters
    ----------
    stmts : list[Statement]
        A list of INDRA statements.
    """

    def pick_ns(ag):
        for ns in ['HGNC', 'FPLX', 'CHEMBL', 'CHEBI', 'GO', 'MESH']:
            if ns in ag.db_refs.keys():
                dbid = ag.db_refs[ns]
                break
        else:
            ns = 'TEXT'
            dbid = ag.name
        return '%s@%s' % (dbid, ns)

    queries = []
    url_base = get_url_base('statements/from_agents')
    non_binary_statements = (Complex, SelfModification, ActiveForm,
                             Translocation, Conversion)
    for stmt in stmts:
        kwargs = {}
        if not isinstance(stmt, non_binary_statements):
            for pos, ag in zip(['subject', 'object'], stmt.agent_list()):
                if ag is not None:
                    kwargs[pos] = pick_ns(ag)
        else:
            for i, ag in enumerate(stmt.agent_list()):
                if ag is not None:
                    kwargs['agent%d' % i] = pick_ns(ag)
        kwargs['type'] = stmt.__class__.__name__
        kwargs.update(params)
        query_str = '?' + '&'.join(['%s=%s' % (k, v) for k, v in kwargs.items()
                                    if v is not None])
        queries.append(url_base + query_str)
    return queries

