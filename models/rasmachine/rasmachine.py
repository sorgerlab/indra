from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import io
import sys
import time
import json
import shutil
import logging
import argparse
import gmail_client
import twitter_client
import ndex.client
from indra import reach
from indra.literature import pubmed_client, get_full_text
from indra.assemblers import CxAssembler, PysbAssembler
from indra.tools.incremental_model import IncrementalModel

model_path = os.path.dirname(os.path.abspath(__file__))
global_filters = ['grounding', 'prior_one', 'human_only']

logger = logging.getLogger('rasmachine')

def get_email_pmids(cred_file, num_days=10):
    try:
        with open(cred_file, 'rt') as fh:
            uname, passwd = [l.strip() for l in fh.readlines()]
    except IOError:
        logger.error('Could not access Gmail credentials at %s.' % cred_file)
        return []

    M = gmail_client.gmail_login(uname, passwd)
    gmail_client.select_mailbox(M, 'INBOX')
    pmids = gmail_client.get_message_pmids(M, num_days)
    return pmids

def get_searchterm_pmids(search_terms, num_days=1):
    pmids = set([])
    for s in search_terms:
        ids = pubmed_client.get_ids(s, reldate=num_days)
        pmids = pmids.union(ids)
    return list(pmids)

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
        if txt_format == 'nxml':
            rp = reach.process_nxml_str(txt, citation=pmid, offline=True)
            if os.path.exists('reach_output.json'):
                shutil.move('reach_output.json', fulltext_path)
        elif txt_format == 'txt':
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
    for pmid in pmids:
        # If the paper has not been included in the model yet
        if model.stmts.get(pmid) is None:
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

def upload_to_ndex(stmts, cred_file):
    try:
        with open(cred_file, 'rt') as fh:
            uname, passwd, network_id = [l.strip() for l in fh.readlines()]
    except IOError:
        logger.error('Could not access NDEx credentials.')
        return
    nd = ndex.client.Ndex('http://public.ndexbio.org',
                            username=uname, password=passwd)
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', help='Model name', required=True)
    parser.add_argument('--twitter', help='Twitter credentials file')
    parser.add_argument('--gmail', help='Gmail credentials file')
    parser.add_argument('--ndex', help='NDEx credentials file')
    parser.add_argument('--belief', help='Belief threshold (between 0 and 1')
    args = parser.parse_args()

    logger.info('-------------------------')
    logger.info(time.strftime('%c'))

    if not args.model:
        logger.error('Model name must be supplied as --model model_name.')
        sys.exit()
    else:
        model_name = args.model

    if args.twitter:
        twitter_cred = args.twitter
        if os.path.exists(twitter_cred):
            use_twitter = True
            logger.info('Using Twitter with credentials: %s' % twitter_cred)
        else:
            use_twitter = False
    else:
        use_twitter = False

    if args.gmail:
        gmail_cred = args.gmail
        if os.path.exists(gmail_cred):
            use_gmail = True
            logger.info('Using Gmail with credentials: %s' % gmail_cred)
        else:
            use_gmail = False
    else:
        use_gmail = False

    if args.ndex:
        ndex_cred = args.ndex
        if os.path.exists(ndex_cred):
            use_ndex = True
            logger.info('Using NDEx with credentials: %s' % ndex_cred)
        else:
            use_ndex = False
    else:
        use_ndex = False

    # Probability cutoff for filtering statements
    if args.belief:
        if not os.path.exists(args.belief):
            BELIEF_THRESHOLD = 0.95
        with open(args.belief, 'rt') as f:
            belief_str = f.read().strip()
        BELIEF_THRESHOLD = float(belief_str)
    else:
        BELIEF_THRESHOLD = 0.95
    logger.info('Using belief threshold: %.2f' % BELIEF_THRESHOLD)

    pmids = []
    # Get email PMIDs
    if use_gmail:
        logger.info('Getting PMIDs from emails.')
        try:
            email_pmids = get_email_pmids(gmail_cred, num_days=10)
            logger.info('Collected %d PMIDs from Gmail' % len(email_pmids))
            pmids += email_pmids
        except Exception as e:
            logger.error('Could not get PMIDs from Gmail, continuing.')
            logger.error(e)

    # Get search PMIDs
    search_terms_file = os.path.join(model_path, model_name, 'search_terms.txt')
    if os.path.exists(search_terms_file):
        with open(search_terms_file, 'rt') as f:
            search_terms = [l.strip() for l in f.readlines()]
        if search_terms:
            logger.info('Using search terms: %s' % ', '.join(search_terms))
            search_pmids = get_searchterm_pmids(search_terms, num_days=5)
            logger.info('Collected %d PMIDs from PubMed' % len(search_pmids))
            pmids += search_pmids

    if not pmids:
        logger.info('No PMIDs found.')
        sys.exit()
    else:
        logger.info('%s PMIDs found in total.' % len(pmids))

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
        elif s.belief > BELIEF_THRESHOLD:
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
        elif s.belief > BELIEF_THRESHOLD:
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
