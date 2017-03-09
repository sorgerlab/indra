from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import io
import sys
import yaml
import time
import json
import shutil
import logging
import datetime
import argparse
import gmail_client
import twitter_client
import ndex.client
from indra import reach
import indra.tools.assemble_corpus as ac
from indra.tools.gene_network import GeneNetwork
from indra.literature import pubmed_client, get_full_text, elsevier_client
from indra.assemblers import CxAssembler, PysbAssembler
from indra.tools.incremental_model import IncrementalModel

model_path = os.path.dirname(os.path.abspath(__file__))
global_filters = ['grounding', 'prior_one', 'human_only']

logger = logging.getLogger('rasmachine')

def build_prior(genes, out_file):
    gn = GeneNetwork(genes)
    stmts = gn.get_statements(filter=False)
    ac.dump_statements(stmts, out_file)
    return stmts

def get_email_pmids(gmail_cred, num_days=10):
    M = gmail_client.gmail_login(gmail_cred.get('user'),
                                 gmail_cred.get('password'))
    gmail_client.select_mailbox(M, 'INBOX')
    pmids = gmail_client.get_message_pmids(M, num_days)
    return pmids

def get_searchterm_pmids(search_terms, num_days=1):
    pmids = {}
    for s in search_terms:
        # Special cases
        if s.upper() == 'MET':
            s = 'c-MET'
        elif s.upper() == 'JUN':
            s = 'c-JUN'
        ids = pubmed_client.get_ids(s, reldate=num_days)
        pmids[s] = ids
    return pmids

def get_searchgenes_pmids(search_genes, num_days=1):
    pmids = {}
    for s in search_genes:
        try:
            ids = pubmed_client.get_ids_for_gene(s, reldate=num_days)
        except ValueError as e:
            logger.error('Gene symbol %s is invalid')
        pmids[s] = ids
    return pmids

def check_pmids(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid is not None:
                if not ev.pmid.isdigit():
                    logger.warning('Invalid PMID: %s' % ev.pmid)

def process_paper(model_name, pmid):
    json_path = os.path.join(model_path, model_name,
                             'jsons', 'PMID%s.json' % pmid)

    if pmid.startswith('api') or pmid.startswith('PMID'):
        logger.warning('Invalid PMID: %s' % pmid)
    # If the paper has been read, use the json output file
    if os.path.exists(json_path):
        rp = reach.process_json_file(json_path, citation=pmid)
        txt_format = 'existing_json'
    # If the paper has not been read, download the text and read
    else:
        txt, txt_format = get_full_text(pmid, 'pmid')
        if txt_format == 'pmc_oa_xml':
            rp = reach.process_nxml_str(txt, citation=pmid, offline=True)
            if os.path.exists('reach_output.json'):
                shutil.move('reach_output.json', json_path)
        elif txt_format == 'elsevier_xml':
            # Extract the raw text from the Elsevier XML
            txt = elsevier_client.extract_text(txt)
            rp = reach.process_text(txt, citation=pmid, offline=True)
            if os.path.exists('reach_output.json'):
                shutil.move('reach_output.json', json_path)
        elif txt_format == 'abstract':
            rp = reach.process_text(txt, citation=pmid, offline=True)
            if os.path.exists('reach_output.json'):
                shutil.move('reach_output.json', json_path)
        else:
            rp = None
    if rp is not None:
        check_pmids(rp.statements)
    return rp, txt_format

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

def extend_model(model_name, model, pmids):
    npapers = 0
    nabstracts = 0
    nexisting = 0
    id_used = []
    for search_term, pmid_list in pmids.items():
        for pmid in pmid_list:
            # If the paper has not been included in the model yet
            if pmid in id_used:
                continue
            id_used.append(pmid)
            if model.stmts.get(pmid) is None:
                logger.info('Processing %s for search term %s' % \
                            (pmid, search_term))
                rp, txt_format = process_paper(model_name, pmid)
                if rp is not None:
                    if txt_format == 'abstract':
                        nabstracts += 1
                    elif txt_format in ['pmc_oa_xml', 'elsevier_xml']:
                        npapers += 1
                    else:
                        nexisting += 1
                    if not rp.statements:
                        logger.info('No statement from PMID%s (%s)' % \
                                    (pmid, txt_format))
                    else:
                        logger.info('%d statements from PMID%s (%s)' % \
                                    (len(rp.statements), pmid, txt_format))
                    model.add_statements(pmid, rp.statements)
                else:
                    logger.info('Reach processing failed for PMID%s' % pmid)
    return npapers, nabstracts, nexisting

def _increment_ndex_ver(ver_str):
    if not ver_str:
        new_ver = '1.0'
    else:
        major_ver, minor_ver = ver_str.split('.')
        new_minor_ver = str(int(minor_ver) + 1)
        new_ver = major_ver + '.' + new_minor_ver
    return new_ver

def upload_to_ndex(stmts, ndex_cred):
    nd = ndex.client.Ndex('http://public.ndexbio.org',
                          username=ndex_cred.get('user'),
                          password=ndex_cred.get('password'))
    network_id = ndex_cred.get('network')

    ca = CxAssembler()
    ca.network_name = 'rasmachine'
    ca.add_statements(stmts)
    ca.make_model()
    cx_str = ca.print_cx()

    try:
        logger.info('Getting network summary...')
        summary = nd.get_network_summary(network_id)
    except Exception as e:
        logger.error('Could not get NDEx network summary.')
        logger.error(e)
        return

    try:
        logger.info('Updating network...')
        cx_stream = io.BytesIO(cx_str.encode('utf-8'))
        with open('cx_upload.cx', 'wb') as fh:
            fh.write(cx_str.encode('utf-8'))
        nd.update_cx_network(cx_stream, network_id)
    except Exception as e:
        logger.error('Could not update NDEx network.')
        logger.error(e)
        return
    ver_str = summary.get('version')
    new_ver = _increment_ndex_ver(ver_str)
    profile = {'name': summary.get('name'),
               'description': summary.get('description'),
               'version': new_ver,
               }
    logger.info('Updating NDEx network (%s) profile to %s' %
                (network_id, profile))
    profile_retries = 5
    for i in range(profile_retries):
        try:
            time.sleep(5)
            nd.update_network_profile(network_id, profile)
            break
        except Exception as e:
            logger.error('Could not update NDEx network profile.')
            logger.error(e)

class InvalidConfigurationError(Exception):
    pass

def get_config(config_fname):
    try:
        fh = open(config_fname, 'rt')
    except IOError as e:
        logger.error('Could not open configuration file %s.' % config_fname)
        raise(e)
    try:
        config = yaml.load(fh)
    except Exception as e:
        logger.error('Could not parse YAML configuration %s.' % config_name)
        raise(e)

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
            if st.belief >= belief_cutoff or \
                sources.intersection(db_names):
                supp.append(st)
        for st in stmt.supported_by:
            sources = set([ev.source_api for ev in st.evidence])
            if st.belief >= belief_cutoff or \
                sources.intersection(db_names):
                supp_by.append(st)
        stmt.supports = supp
        stmt.supported_by = supp_by
    logger.info('%d statements after filter...' % len(stmts_out))
    return stmts_out

if __name__ == '__main__':
    logger.info('-------------------------')
    logger.info(time.strftime('%c'))

    if len(sys.argv) < 2:
        logger.error('Model name argument missing')
        sys.exit()
    if len(sys.argv) < 3:
        logger.error('Configuration file argument missing.')
        sys.exit()

    model_name = sys.argv[1]
    config_fname = sys.argv[2]
    try:
        config = get_config(config_fname)
    except Exception as e:
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

    twitter_cred = config.get('twitter')
    if twitter_cred:
        use_twitter = True
        if not twitter_cred.get('consumer_token'):
            logger.info('Twitter consumer token (consumer_token) missing.')
            use_twitter = False
        if not twitter_cred.get('consumer_secret'):
            logger.info('Twitter consumer secret (consumer_secret) missing.')
            use_twitter = False
        if not twitter_cred.get('access_token'):
            logger.info('Twitter access token (access_token) missing.')
            use_twitter = False
        if not twitter_cred.get('access_secret'):
            logger.info('Twitter access secret (access_secret) missing.')
            use_twitter = False
    else:
        use_twitter = False
    if use_twitter:
        logger.info('Using Twitter with given credentials.')
    else:
        logger.info('Not using Twitter due to missing credentials.')

    gmail_cred = config.get('gmail')
    if gmail_cred:
        use_gmail = True
        if not gmail_cred.get('user'):
            logger.info('Gmail user missing.')
            use_gmail = False
        if not gmail_cred.get('password'):
            logger.info('Gmail password missing.')
            use_gmail = False
    else:
        use_gmail = False
    if use_gmail:
        logger.info('Using Gmail with given credentials.')
    else:
        logger.info('Not using Gmail due to missing credentials.')

    ndex_cred = config.get('ndex')
    if ndex_cred:
        use_ndex = True
        if not ndex_cred.get('user'):
            logger.info('NDEx user missing.')
            use_ndex = False
        if not ndex_cred.get('password'):
            logger.info('NDEx password missing.')
            use_ndex = False
        if not ndex_cred.get('network'):
            logger.info('NDEx network missing.')
            use_ndex = False
    if use_ndex:
        logger.info('Using NDEx with given credentials.')
    else:
        logger.info('Not using NDEx due to missing information.')

    pmids = {}
    # Get email PMIDs
    if use_gmail:
        logger.info('Getting PMIDs from emails.')
        try:
            email_pmids = get_email_pmids(gmail_cred, num_days=10)
            # Put the email_pmids into the pmids dictionary
            pmids['Gmail'] = email_pmids
            logger.info('Collected %d PMIDs from Gmail' % len(email_pmids))
        except Exception as e:
            logger.error('Could not get PMIDs from Gmail, continuing.')
            logger.error(e)

    # Get PMIDs for general search_terms
    search_terms = config.get('search_terms')
    if not search_terms:
        logger.info('No search terms argument (search_terms) specified.')
    else:
        search_genes = config.get('search_genes')
        if search_genes:
            search_terms += search_genes
        logger.info('Using search terms: %s' % ', '.join(search_terms))
        pmids_term = get_searchterm_pmids(search_terms, num_days=5)
        num_pmids = sum([len(pm) for pm in pmids_term.values()])
        logger.info('Collected %d PMIDs from PubMed search_terms.' % num_pmids)
        pmids = _extend_dict(pmids, pmids_term)


    search_genes = config.get('search_genes')

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


    # Load the model
    logger.info(time.strftime('%c'))
    logger.info('Loading original model.')
    inc_model_file = os.path.join(model_path, model_name, 'model.pkl')
    model = IncrementalModel(inc_model_file)
    # Include search genes as prior genes
    model.prior_genes = search_genes
    stats = {}
    logger.info(time.strftime('%c'))
    logger.info('Preassembling original model.')
    model.preassemble(filters=global_filters)
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
                            extend_model(model_name, model, pmids)
    # Having added new statements, we preassemble the model
    model.preassemble(filters=global_filters)

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
    date = datetime.datetime.today()
    date_str = date.strftime('%Y-%m-%d-%H-%M-%S')
    inc_model_bkp_file = os.path.join(model_path, model_name,
                                      'model-%s.pkl' % date_str)
    model.save(inc_model_bkp_file)

    # Upload the new, final statements to NDEx
    if use_ndex:
        logger.info('Uploading to NDEx')
        logger.info(time.strftime('%c'))
        upload_to_ndex(new_stmts, ndex_cred)

    # Print and tweet the status message
    logger.info('--- Final statistics ---')
    for k, v in sorted(stats.items(), key=lambda x: x[0]):
        logger.info('%s: %s' % (k, v))
    logger.info('------------------------')

    msg_str = make_status_message(stats)
    if msg_str is not None:
        logger.info('Status message: %s' % msg_str)
        if use_twitter:
            logger.info('Now tweeting: %s' % msg_str)
            twitter_client.update_status(msg_str, twitter_cred)
