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
import argparse
import gmail_client
import twitter_client
import ndex.client
from indra import reach
from indra.literature import pubmed_client, get_full_text, elsevier_client
from indra.assemblers import CxAssembler, PysbAssembler
from indra.tools.incremental_model import IncrementalModel

model_path = os.path.dirname(os.path.abspath(__file__))
global_filters = ['grounding', 'prior_one', 'human_only']

logger = logging.getLogger('rasmachine')

def get_email_pmids(gmail_cred, num_days=10):
    M = gmail_client.gmail_login(gmail_cred.get('user'),
                                 gmail_cred.get('password'))
    gmail_client.select_mailbox(M, 'INBOX')
    pmids = gmail_client.get_message_pmids(M, num_days)
    return pmids

def get_searchterm_pmids(search_terms, num_days=1):
    pmids = {}
    for s in search_terms:
        ids = pubmed_client.get_ids(s, reldate=num_days)
        pmids[s] = ids
    return pmids

def process_paper(model_name, pmid):
    abstract_path = os.path.join(model_path, model_name,
                                 'jsons', 'abstract', 'PMID%s.json' % pmid)
    fulltext_path = os.path.join(model_path, model_name,
                                 'jsons', 'full', 'PMID%s.json' % pmid)

    # If the paper has been parsed, use the parse output file
    if os.path.exists(abstract_path):
        rp = reach.process_json_file(abstract_path, citation=pmid)
        txt_format = 'abstract'
    elif os.path.exists(fulltext_path):
        rp = reach.process_json_file(fulltext_path, citation=pmid)
        txt_format = 'txt'
    # If the paper has not been parsed, download the text and parse
    else:
        txt, txt_format = get_full_text(pmid, 'pmid')
        if txt_format == 'pmc_oa_xml':
            rp = reach.process_nxml_str(txt, citation=pmid, offline=True)
            if os.path.exists('reach_output.json'):
                shutil.move('reach_output.json', fulltext_path)
        elif txt_format == 'elsevier_xml':
            # Extract the raw text from the Elsevier XML
            txt = elsevier_client.extract_text(txt)
            rp = reach.process_text(txt, citation=pmid, offline=True)
            if os.path.exists('reach_output.json'):
                shutil.move('reach_output.json', fulltext_path)
        elif txt_format == 'abstract':
            rp = reach.process_text(txt, citation=pmid, offline=True)
            if os.path.exists('reach_output.json'):
                shutil.move('reach_output.json', abstract_path)
        else:
            rp = None
    return rp, txt_format

def make_status_message(stats):
    ndiff = (stats['new_top'] - stats['orig_top'])
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
                    else:
                        npapers += 1
                    if not rp.statements:
                        logger.info('No statement from PMID%s (%s)' % \
                                    (pmid, txt_format))
                    else:
                        logger.info('%d statements from PMID%s (%s)' % \
                                    (len(rp.statements), pmid, txt_format))
                    model.add_statements(pmid, rp.statements)
                else:
                    logger.info('Reach processing failed for PMID%s' % pmid)
    # Having added new statements, we preassemble the model
    # to merge duplicated and find related statements
    model.preassemble(filters=global_filters)
    return npapers, nabstracts

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
        summary = nd.get_network_summary(network_id)
    except Exception as e:
        logger.error('Could not get NDEx network summary.')
        logger.error(e)
        return

    try:
        nd.update_cx_network(cx_str, network_id)
    except Exception as e:
        logger.error('Could not update NDEx network.')
        logger.error(e)
        return
    ver_str = summary.get('version')
    new_ver = _increment_ndex_ver(ver_str)
    profile = {'name': summary.get('name'),
               'description': summary.get('description'),
               'version': new_ver,
               'visibility': 'PUBLIC'
               }
    try:
        nd.update_network_profile(network_id, profile)
    except Exception as e:
        logger.error('Could not update NDEx network profile.')
        logger.error(e)
        return

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

    email_pmids = []
    # Get email PMIDs
    if use_gmail:
        logger.info('Getting PMIDs from emails.')
        try:
            email_pmids = get_email_pmids(gmail_cred, num_days=10)
            logger.info('Collected %d PMIDs from Gmail' % len(email_pmids))
        except Exception as e:
            logger.error('Could not get PMIDs from Gmail, continuing.')
            logger.error(e)

    pmids = {}
    search_terms = config.get('search_terms')
    if not search_terms:
        logging.info('No search terms argument (search_terms) specified.')
    else:
        logger.info('Using search terms: %s' % ', '.join(search_terms))
        pmids = get_searchterm_pmids(search_terms, num_days=5)
        num_pmids = sum([len(pm) for _, pm in pmids.items()])
        logger.info('Collected %d PMIDs from PubMed.' % num_pmids)

    # Put the email_pmids into the pmids dictionary
    pmids['Gmail'] = email_pmids

    # Load the model
    logger.info(time.strftime('%c'))
    logger.info('Loading original model.')
    inc_model_file = os.path.join(model_path, model_name, 'model.pkl')
    model = IncrementalModel(inc_model_file)
    stats = {}
    logger.info(time.strftime('%c'))
    logger.info('Preassembling original model.')
    model.preassemble(filters=global_filters)
    logger.info(time.strftime('%c'))

    # Original statistics
    stats['orig_stmts'] = len(model.get_statements())
    stats['orig_unique'] = len(model.unique_stmts)
    stats['orig_top'] = len(model.toplevel_stmts)
    # Filter the top level statements with a probability cutoff
    orig_likely = []
    for s in model.toplevel_stmts:
        if s.evidence[0].source_api in ['bel', 'biopax']:
            orig_likely.append(s)
        elif s.belief > belief_threshold:
            orig_likely.append(s)
    stats['orig_likely'] = len(orig_likely)

    # Make a PySB model from filtered statements
    #pysb_assmb = PysbAssembler()
    #pysb_assmb.add_statements(orig_likely)
    #pysb_assmb.make_model()
    # Stats for Pysb assembled model
    #stats['orig_monomers'] = len(pysb_assmb.model.monomers)
    #stats['orig_rules'] = len(pysb_assmb.model.rules)

    # Extend the model with PMIDs
    logger.info('----------------')
    logger.info(time.strftime('%c'))
    logger.info('Extending model.')
    stats['new_papers'], stats['new_abstracts'] = \
                            extend_model(model_name, model, pmids)

    # New statistics
    stats['new_stmts'] = len(model.get_statements())
    stats['new_unique'] = len(model.unique_stmts)
    stats['new_top'] = len(model.toplevel_stmts)
    new_likely = []
    for s in model.toplevel_stmts:
        if s.evidence[0].source_api in ['bel', 'biopax']:
            new_likely.append(s)
        elif s.belief > belief_threshold:
            new_likely.append(s)
    stats['new_likely'] = len(new_likely)

    # Make a PySB model from filtered statements
    #pysb_assmb = PysbAssembler()
    #pysb_assmb.add_statements(new_likely)
    #pysb_assmb.make_model()
    # Stats for Pysb assembled model
    #stats['new_monomers'] = len(pysb_assmb.model.monomers)
    #stats['new_rules'] = len(pysb_assmb.model.rules)

    # Save model
    logger.info(time.strftime('%c'))
    logger.info('Saving model')
    model.save(inc_model_file)
    logger.info(time.strftime('%c'))

    # Upload the new, highly likely statements to NDEx
    if use_ndex:
        logger.info('Uploading to NDEx')
        logger.info(time.strftime('%c'))
        upload_to_ndex(new_likely, ndex_cred)

    # Print and tweet the status message
    logger.info('--- Final statistics ---')
    for k, v in stats.items():
        logger.info('%s: %s' % (k, v))
    logger.info('------------------------')

    msg_str = make_status_message(stats)
    if msg_str is not None:
        logger.info('Status message: %s' % msg_str)
        if use_twitter:
            logger.info('Now tweeting: %s' % msg_str)
            twitter_client.update_status(msg_str, twitter_cred)
