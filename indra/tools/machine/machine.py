from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import sys
import yaml
import time
import click
import shutil
import tzlocal
import logging
import datetime
import itertools as itt
from indra.sources import reach
from collections import defaultdict
from indra.databases import ndex_client
import indra.tools.assemble_corpus as ac
from indra.assemblers import CxAssembler
from indra.tools.machine import gmail_client
from indra.tools.machine import twitter_client
from indra.tools.gene_network import GeneNetwork
from indra.tools.incremental_model import IncrementalModel
from indra.literature import pubmed_client, get_full_text, elsevier_client
from indra import get_config, has_config

try:
    import boto3
    from indra.tools.reading.submit_reading_pipeline import \
        submit_reading, wait_for_complete
    # Try to make a client
    client = boto3.client('batch')
    from indra.literature.s3_client import get_reader_json_str, get_full_text
    aws_available = True
except Exception:
    aws_available = False

global_filters = ['grounding', 'prior_one', 'human_only']

logger = logging.getLogger('rasmachine')


def build_prior(genes, out_file):
    gn = GeneNetwork(genes)
    stmts = gn.get_statements(filter=False)
    ac.dump_statements(stmts, out_file)
    return stmts


def get_email_pmids(gmail_cred):
    mailbox = gmail_client.gmail_login(gmail_cred.get('user'),
                                 gmail_cred.get('password'))
    gmail_client.select_mailbox(mailbox, 'INBOX')
    num_days = int(gmail_cred.get('num_days', 10))
    logger.info('Searching last %d days of emails', num_days)
    pmids = gmail_client.get_message_pmids(mailbox, num_days)
    return pmids


def get_searchterm_pmids(search_terms, num_days):
    pmids = {}
    for s in search_terms:
        # Special cases
        if s.upper() == 'MET':
            s = 'c-MET'
        elif s.upper() == 'JUN':
            s = 'c-JUN'
        pmids[s] = pubmed_client.get_ids(s, reldate=num_days)
    return pmids


def get_searchgenes_pmids(search_genes, num_days):
    pmids = {}
    for s in search_genes:
        try:
            pmids[s] = pubmed_client.get_ids_for_gene(s, reldate=num_days)
        except ValueError:
            logger.error('Gene symbol %s is invalid')
            continue
    return pmids


def check_pmids(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid is not None:
                if not ev.pmid.isdigit():
                    logger.warning('Invalid PMID: %s' % ev.pmid)


def process_paper(model_name, pmid):
    """Process a paper with the given pubmed identifier

    Parameters
    ----------
    model_name : str
        The directory for the INDRA machine
    pmid : str
        The PMID to process.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    txt_format : str
        A string representing the format of the text
    """
    json_directory = os.path.join(model_name, 'jsons')
    json_path = os.path.join(json_directory, 'PMID%s.json' % pmid)

    if pmid.startswith('api') or pmid.startswith('PMID'):
        logger.warning('Invalid PMID: %s' % pmid)
    # If the paper has been read, use the json output file
    if os.path.exists(json_path):
        rp = reach.process_json_file(json_path, citation=pmid)
        txt_format = 'existing_json'
    # If the paper has not been read, download the text and read
    else:
        try:
            txt, txt_format = get_full_text(pmid, 'pmid')
        except Exception:
            return None, None

        if txt_format == 'pmc_oa_xml':
            rp = reach.process_nxml_str(txt, citation=pmid, offline=True,
                                        output_fname=json_path)
        elif txt_format == 'elsevier_xml':
            # Extract the raw text from the Elsevier XML
            txt = elsevier_client.extract_text(txt)
            rp = reach.process_text(txt, citation=pmid, offline=True,
                                    output_fname=json_path)
        elif txt_format == 'abstract':
            rp = reach.process_text(txt, citation=pmid, offline=True,
                                    output_fname=json_path)
        else:
            rp = None
    if rp is not None:
        check_pmids(rp.statements)
    return rp, txt_format


def process_paper_aws(pmid, start_time_local):
    try:
        metadata, content_type = get_full_text(pmid, metadata=True)
    except Exception as e:
        logger.error('Could not get content from S3: %s' % e)
        return None, None
    logger.info('Downloading %s output from AWS' % pmid)
    reach_json_str = get_reader_json_str('reach', pmid)
    if not reach_json_str:
        logger.info('Could not get output.')
        return None, content_type
    rp = reach.process_json_str(reach_json_str)

    current_time_local = datetime.datetime.now(tzlocal.get_localzone())
    dt_script = current_time_local - start_time_local
    last_mod_remote = metadata['LastModified']
    dt = (current_time_local - last_mod_remote)
    # If it was not modified since the script started
    if dt > dt_script:
        content_type = 'existing_json'
    return rp, content_type


def make_status_message(stats):
    ndiff = (stats['new_final'] - stats['orig_final'])
    msg_str = None
    if (((stats['new_papers'] > 0) or
        (stats['new_abstracts'] > 0)) and
        (ndiff > 0)):
        papers_str = '%d paper' % stats['new_papers']
        if stats['new_papers'] > 1:
            papers_str += 's'
        abstr_str = '%d abstract' % stats['new_abstracts']
        if stats['new_abstracts'] > 1:
            abstr_str += 's'
        mech_str = '%d new mechanism' % ndiff
        if ndiff > 1:
            mech_str += 's'

        if stats['new_papers'] > 0:
            if stats['new_abstracts'] > 0:
                msg_str = 'Today I read %s and %s, and learned %s!' %\
                    (papers_str, abstr_str, mech_str)
            else:
                msg_str = 'Today I read %s, and learned %s!' %\
                    (papers_str, mech_str)
        else:
            if stats['new_abstracts'] > 0:
                msg_str = 'Today I read %s, and learned %s!' %\
                    (abstr_str, mech_str)
    return msg_str


def process_paper_helper(model_name, pmid, start_time_local):
    """Wraps processing a paper by either a local or remote service
    and caches any uncaught exceptions"""
    try:
        if not aws_available:
            rp, txt_format = process_paper(model_name, pmid)
        else:
            rp, txt_format = process_paper_aws(pmid, start_time_local)
    except:
        logger.exception('uncaught exception while processing %s', pmid)
        return None, None

    return rp, txt_format


def extend_model(model_name, model, pmids, start_time_local):
    npapers = 0
    nabstracts = 0
    nexisting = 0

    # Preprocess PMID search results
    pmids_inv = defaultdict(list)
    for search_term, pmid_list in pmids.items():
        for pmid in pmid_list:
            if pmid in model.stmts:
                continue
            pmids_inv[pmid].append(search_term)

    logger.info('Found %d unique and novel PMIDS', len(pmids_inv))

    if not os.path.exists(os.path.join(model_name, 'jsons')):
        os.mkdir(os.path.join(model_name, 'jsons'))

    for counter, (pmid, search_terms) in enumerate(pmids_inv.items(), start=1):
        logger.info('[%d/%d] Processing %s for search terms: %s',
                    counter, len(pmids_inv), pmid, search_terms)

        rp, txt_format = process_paper_helper(model_name, pmid,
                                              start_time_local)

        if rp is None:
            logger.info('Reach processing failed for PMID%s', pmid)
            continue

        if txt_format == 'abstract':
            nabstracts += 1
        elif txt_format in ['pmc_oa_xml', 'elsevier_xml']:
            npapers += 1
        else:
            nexisting += 1

        if not rp.statements:
            logger.info('No statement from PMID%s (%s)' %
                        (pmid, txt_format))
        else:
            logger.info('%d statements from PMID%s (%s)' %
                        (len(rp.statements), pmid, txt_format))
        model.add_statements(pmid, rp.statements)

    return npapers, nabstracts, nexisting


def assemble_cx(stmts, name):
    ca = CxAssembler()
    ca.network_name = name
    ca.add_statements(stmts)
    ca.make_model()
    cx_str = ca.print_cx()
    return cx_str


class InvalidConfigurationError(Exception):
    pass


def get_machine_config(config_fname):
    try:
        fh = open(config_fname, 'rt')
    except IOError as e:
        logger.error('Could not open configuration file %s.' % config_fname)
        raise e
    try:
        config = yaml.load(fh)
    except Exception as e:
        logger.error('Could not parse YAML configuration %s.' % config_fname)
        raise e

    return config


def _extend_dict(d1, d2):
    for k, v in d2.items():
        d1[k] = v
    return d1


def filter_db_highbelief(stmts_in, db_names, belief_cutoff):
    logger.info('Filtering %d statements to above %f belief' %
                (len(stmts_in), belief_cutoff))
    # The first round of filtering is in the top-level list
    stmts_out = []
    # Now we eliminate supports/supported-by
    for stmt in stmts_in:
        sources = set([ev.source_api for ev in stmt.evidence])
        if stmt.belief >= belief_cutoff or \
            sources.intersection(db_names):
            stmts_out.append(stmt)
        else:
            continue
        supp_by = []
        supp = []
        for st in stmt.supports:
            sources = set([ev.source_api for ev in st.evidence])
            if (st.belief >= belief_cutoff or
                sources.intersection(db_names)):
                supp.append(st)
        for st in stmt.supported_by:
            sources = set([ev.source_api for ev in st.evidence])
            if (st.belief >= belief_cutoff or
               sources.intersection(db_names)):
                supp_by.append(st)
        stmt.supports = supp
        stmt.supported_by = supp_by
    logger.info('%d statements after filter...' % len(stmts_out))
    return stmts_out


def upload_new_ndex(model_path, new_stmts, ndex_cred):
    logger.info('Uploading to NDEx')
    logger.info(time.strftime('%c'))
    cx_str = assemble_cx(new_stmts, name=ndex_cred.get('name', 'rasmachine'))
    cx_name = os.path.join(model_path, 'model.cx')
    with open(cx_name, 'wb') as fh:
        fh.write(cx_str.encode('utf-8'))
    upload_cx_to_ndex(cx_str, ndex_cred)


def upload_cx_to_ndex(cx_str, ndex_cred):
    network_id = ndex_cred['network']
    provenance = make_ndex_provenance(network_id)
    ndex_client.set_provenance(provenance, network_id, ndex_cred)
    try:
        ndex_client.update_network(cx_str, network_id, ndex_cred)
    except Exception as e:
        logger.error('NDEx network update failed')
        logger.exception(e)


def make_ndex_provenance(network_id):
    ts = int(time.time())
    uri = 'http://public.ndexbio.org/v2/network/%s/summary' % network_id
    inputs = [{'creationEvent':
                {'eventType': 'Create Network',
                 'inputs': []},
               'uri': uri,
               'properties': []
              }]
    properties = [{
        'name': 'update_count',
        'value': '400'
        }]
    creation_event = {
        'eventType': 'Update Network',
        'eventAtTime': ts,
        'startedAtTime': 1464633540,
        'inputs': inputs,
        'uri': uri,
        'properties': properties
        }
    return creation_event


def make_date_str():
    return datetime.datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S')


def get_gmail_cred(config):
    gmail_cred = config.get('gmail')
    if not gmail_cred:
        return
    elif not gmail_cred.get('user'):
        logger.info('Gmail user missing.')
        return
    elif not gmail_cred.get('password'):
        logger.info('Gmail password missing.')
        return
    return gmail_cred


def get_twitter_cred(config):
    twitter_cred = config.get('twitter')

    if not twitter_cred:
        return
    elif not twitter_cred.get('consumer_token'):
        logger.info('Twitter consumer token (consumer_token) missing.')
        return
    elif not twitter_cred.get('consumer_secret'):
        logger.info('Twitter consumer secret (consumer_secret) missing.')
        return
    elif not twitter_cred.get('access_token'):
        logger.info('Twitter access token (access_token) missing.')
        return
    elif not twitter_cred.get('access_secret'):
        logger.info('Twitter access secret (access_secret) missing.')
        return

    return twitter_cred


def get_ndex_cred(config):
    ndex_cred = config.get('ndex')
    if not ndex_cred:
        return
    elif not ndex_cred.get('user') and not has_config('NDEX_USERNAME'):
        logger.info('NDEx user missing.')
        return
    elif not ndex_cred.get('password')and not has_config('NDEX_PASSWORD'):
        logger.info('NDEx password missing.')
        return
    elif not ndex_cred.get('network'):
        logger.info('NDEx network missing.')
        return
    return ndex_cred


def run_machine(model_path, pmids, belief_threshold, search_genes=None,
                ndex_cred=None, twitter_cred=None, grounding_map=None):
    start_time_local = datetime.datetime.now(tzlocal.get_localzone())
    date_str = make_date_str()

    # Save PMIDs in file and send for remote reading
    if aws_available:
        pmid_fname = 'pmids-%s.txt' % date_str
        all_pmids = []
        for v in pmids.values():
            all_pmids += v
        all_pmids = list(set(all_pmids))

        with open(pmid_fname, 'wt') as fh:
            for pmid in all_pmids:
                fh.write('%s\n' % pmid)
        # Submit reading
        job_list = submit_reading('rasmachine', pmid_fname, ['reach'])

        # Wait for reading to complete
        wait_for_complete('run_reach_queue', job_list, idle_log_timeout=600,
                          kill_on_log_timeout=True)

    # Load the model
    logger.info(time.strftime('%c'))
    logger.info('Loading original model.')
    inc_model_file = os.path.join(model_path, 'model.pkl')
    model = IncrementalModel(inc_model_file)
    # Include search genes as prior genes
    if search_genes:
        model.prior_genes = search_genes
    stats = {}
    logger.info(time.strftime('%c'))
    logger.info('Preassembling original model.')
    model.preassemble(filters=global_filters, grounding_map=grounding_map)
    logger.info(time.strftime('%c'))

    # Original statistics
    stats['orig_stmts'] = len(model.get_statements())
    stats['orig_assembled'] = len(model.assembled_stmts)
    orig_stmts = filter_db_highbelief(model.assembled_stmts, ['bel', 'biopax'],
                                      belief_threshold)
    orig_stmts = ac.filter_top_level(orig_stmts)
    stats['orig_final'] = len(orig_stmts)
    logger.info('%d final statements' % len(orig_stmts))

    # Extend the model with PMIDs
    logger.info('----------------')
    logger.info(time.strftime('%c'))
    logger.info('Extending model.')
    stats['new_papers'], stats['new_abstracts'], stats['existing'] = \
        extend_model(model_path, model, pmids, start_time_local)
    # Having added new statements, we preassemble the model
    model.preassemble(filters=global_filters, grounding_map=grounding_map)

    # New statistics
    stats['new_stmts'] = len(model.get_statements())
    stats['new_assembled'] = len(model.assembled_stmts)
    new_stmts = filter_db_highbelief(model.assembled_stmts, ['bel', 'biopax'],
                                     belief_threshold)
    new_stmts = ac.filter_top_level(new_stmts)
    stats['new_final'] = len(new_stmts)
    logger.info('%d final statements' % len(new_stmts))

    check_pmids(model.get_statements())

    # Save model
    logger.info(time.strftime('%c'))
    logger.info('Saving model')
    model.save(inc_model_file)
    logger.info(time.strftime('%c'))

    # Save a time stamped version of the pickle for backup/diagnostic purposes
    if not aws_available:
        inc_model_bkp_file = os.path.join(model_path,
                                          'model-%s.pkl' % date_str)
        model.save(inc_model_bkp_file)
    else:
        key = 'rasmachine/%s/model-%s.pkl' % (model_path.replace('/', '_'),
                                              date_str)
        s3 = boto3.client('s3')
        s3.upload_file(inc_model_file, 'bigmech', key)

    # Upload the new, final statements to NDEx
    if ndex_cred:
        upload_new_ndex(model_path, new_stmts, ndex_cred)

    # Print and tweet the status message
    logger.info('--- Final statistics ---')
    for k, v in sorted(stats.items(), key=lambda x: x[0]):
        logger.info('%s: %s' % (k, v))
    logger.info('------------------------')

    msg_str = make_status_message(stats)
    if msg_str is not None:
        logger.info('Status message: %s' % msg_str)
        if twitter_cred:
            logger.info('Now tweeting: %s' % msg_str)
            twitter_client.update_status(msg_str, twitter_cred)


def run_with_search_helper(model_path, config, num_days=None):
    logger.info('-------------------------')
    logger.info(time.strftime('%c'))

    if not os.path.isdir(model_path):
        logger.error('%s is not a directory', model_path)
        sys.exit()

    default_config_fname = os.path.join(model_path, 'config.yaml')

    if config:
        config = get_machine_config(config)
    elif os.path.exists(default_config_fname):
        logger.info('Loading default configuration from %s',
                    default_config_fname)
        config = get_machine_config(default_config_fname)
    else:
        logger.error('Configuration file argument missing.')
        sys.exit()

    # Probability cutoff for filtering statements
    default_belief_threshold = 0.95
    belief_threshold = config.get('belief_threshold')
    if belief_threshold is None:
        belief_threshold = default_belief_threshold
        msg = 'Belief threshold argument (belief_threshold) not specified.' + \
              ' Using default belief threshold %.2f' % default_belief_threshold
        logger.info(msg)
    else:
        logger.info('Using belief threshold: %.2f' % belief_threshold)

    twitter_cred = get_twitter_cred(config)
    if twitter_cred:
        logger.info('Using Twitter with given credentials.')
    else:
        logger.info('Not using Twitter due to missing credentials.')

    gmail_cred = get_gmail_cred(config)
    if gmail_cred:
        logger.info('Using Gmail with given credentials.')
    else:
        logger.info('Not using Gmail due to missing credentials.')

    ndex_cred = get_ndex_cred(config)
    if ndex_cred:
        logger.info('Using NDEx with given credentials.')
    else:
        logger.info('Not using NDEx due to missing information.')

    pmids = {}
    # Get email PMIDs
    if gmail_cred:
        logger.info('Getting PMIDs from emails.')
        try:
            email_pmids = get_email_pmids(gmail_cred)
            # Put the email_pmids into the pmids dictionary
            pmids['Gmail'] = email_pmids
            logger.info('Collected %d PMIDs from Gmail', len(email_pmids))
        except Exception:
            logger.exception('Could not get PMIDs from Gmail, continuing.')

    # Get PMIDs for general search_terms and genes
    search_genes = config.get('search_genes')
    search_terms = config.get('search_terms')
    if not search_terms:
        logger.info('No search terms argument (search_terms) specified.')
    else:
        if search_genes is not None:
            search_terms += search_genes
        logger.info('Using search terms: %s' % ', '.join(search_terms))

        if num_days is None:
            num_days = int(config.get('search_terms_num_days', 5))
        logger.info('Searching the last %d days', num_days)

        pmids_term = get_searchterm_pmids(search_terms, num_days=num_days)
        num_pmids = len(set(itt.chain.from_iterable(pmids_term.values())))
        logger.info('Collected %d PMIDs from PubMed search_terms.', num_pmids)
        pmids = _extend_dict(pmids, pmids_term)

    # Get optional grounding map
    gm_path = config.get('grounding_map_path')
    if gm_path:
        try:
            from indra.preassembler.grounding_mapper import load_grounding_map
            grounding_map = load_grounding_map(gm_path)
        except Exception as e:
            logger.error('Could not load grounding map from %s' % gm_path)
            logger.error(e)
            grounding_map = None
    else:
        grounding_map = None

    '''
    # Get PMIDs for search_genes
    # Temporarily removed because Entrez-based article searches
    # are lagging behind and cannot be time-limited
    if not search_genes:
        logger.info('No search genes argument (search_genes) specified.')
    else:
        logger.info('Using search genes: %s' % ', '.join(search_genes))
        pmids_gene = get_searchgenes_pmids(search_genes, num_days=5)
        num_pmids = sum([len(pm) for pm in pmids_gene.values()])
        logger.info('Collected %d PMIDs from PubMed search_genes.' % num_pmids)
        pmids = _extend_dict(pmids, pmids_gene)
    '''
    run_machine(
        model_path,
        pmids,
        belief_threshold,
        search_genes=search_genes,
        ndex_cred=ndex_cred,
        twitter_cred=twitter_cred,
        grounding_map=grounding_map
    )


def load_model(model_path):
    logger.info(time.strftime('%c'))
    logger.info('Loading original model.')
    inc_model_file = os.path.join(model_path, 'model.pkl')
    return IncrementalModel(inc_model_file)


def summarize_helper(model_path):
    model = load_model(model_path)
    stmts = model.get_statements()
    click.echo('Number of statements: {}'.format(len(stmts)))
    agents = model.get_model_agents()
    click.echo('Number of agents: {}'.format(len(agents)))


def run_with_pmids_helper(model_path, pmids):
    default_config_fname = os.path.join(model_path, 'config.yaml')
    config = get_machine_config(default_config_fname)

    belief_threshold = config.get('belief_threshold', 0.95)
    twitter_cred = get_twitter_cred(config)
    ndex_cred = get_ndex_cred(config)

    # Get optional grounding map
    gm_path = config.get('grounding_map_path')
    if gm_path:
        try:
            from indra.preassembler.grounding_mapper import load_grounding_map
            grounding_map = load_grounding_map(gm_path)
        except Exception as e:
            logger.error('Could not load grounding map from %s' % gm_path)
            logger.error(e)
            grounding_map = None
    else:
        grounding_map = None

    run_machine(
        model_path,
        {'enumerated': [pmid.strip() for pmid in pmids]},
        belief_threshold,
        ndex_cred=ndex_cred,
        twitter_cred=twitter_cred,
        grounding_map=grounding_map
    )
