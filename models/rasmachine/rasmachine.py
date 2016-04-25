import os
import sys
import time
import shutil
import argparse
import gmail_client
import twitter_client
from incremental_model import IncrementalModel
from indra import reach
from indra.literature import pubmed_client, get_full_text

model_path = os.path.dirname(os.path.abspath(__file__))

def get_email_pmids(cred_file):    
    try:
        fh = open(cred_file, 'rt')
        uname, passwd = [l.strip() for l in fh.readlines()]
    except IOError:
        print 'Could not access Gmail credentials.'
        return []

    M = gmail_client.gmail_login(uname, passwd)
    gmail_client.select_mailbox(M, 'INBOX')

    pmids = gmail_client.get_message_pmids(M)
    print 'Collected %d PMIDs' % len(pmids)
    return pmids

def get_searchterm_pmids(search_terms, num_days=1):
    pmids = set([])
    for s in search_terms:
        ids = pubmed_client.get_ids(s, reldate=num_days)
        pmids = pmids.union(ids)
    return list(pmids)

def process_paper(pmid):
    global model_name
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
        txt, txt_format = get_full_text(pmid)
        if txt_format == 'nxml':
            rp = reach.process_nxml_str(txt, citation=pmid, offline=True)
            shutil.move('reach_output.json', fulltext_path)
        elif txt_format == 'txt':
            rp = reach.process_text(txt, citation=pmid, offline=True)
            shutil.move('reach_output.json', fulltext_path)
        elif txt_format == 'abstract':
            rp = reach.process_text(txt, citation=pmid, offline=True)
            shutil.move('reach_output.json', abstract_path)
        else:
            rp = None
    return rp, txt_format

def make_status_message(stats):
    ndiff = (stats['new_stmts'] - stats['orig_stmts'])
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

def extend_model(model, pmids):
    global model_name
    npapers = 0
    nabstracts = 0
    for pmid in pmids:
        # If the paper has not been included in the model yet
        if model.stmts.get(pmid) is None:
            rp, txt_format = process_paper(pmid)
            if rp is not None:
                if txt_format == 'abstract':
                    nabstracts += 1
                else:
                    npapers += 1
                print pmid, len(rp.statements)
                model.add_statements(pmid, rp.statements)
            else:
                model.add_statements(pmid, [])
                print 'No statement extracted from PMID%s' % pmid
    return npapers, nabstracts

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', help='Model name', required=True)
    parser.add_argument('--twitter', help='Twitter credentials file')
    parser.add_argument('--gmail', help='Gmail credentials file')
    args = parser.parse_args()

    print time.strftime('%c')

    global model_name

    if not args.model:
        print 'Model name must be supplied as --model model_name.'
        sys.exit()
    else:
        model_name = args.model

    if args.twitter:
        twitter_cred = args.twitter
        if os.path.exists(twitter_cred):
            use_twitter = True
        else:
            use_twitter = False
    else:
        use_twitter = False

    if args.gmail:
        gmail_cred = args.gmail
        if os.path.exists(gmail_cred):
            use_gmail = True
        else:
            use_gmail = False
    else:
        use_gmail = False

    pmids = []
    # Get email PMIDs
    if use_gmail: 
        pmids += get_email_pmids(gmail_cred)

    # Get search PMIDs
    search_terms_file = os.path.join(model_path, model_name, 'search_terms.txt')
    if os.path.exists(search_terms_file):
        search_terms = [l.strip() for l in
                    open(search_terms_file, 'rt').readlines()]
        if search_terms:
            pmids += get_searchterm_pmids(search_terms) 
    if not pmids:
        print 'No PMIDs found.'
        sys.exit()

    # Load the model
    inc_model_file = os.path.join(model_path, model_name, 'model.pkl')
    model = IncrementalModel(inc_model_file)
    pysb_model = model.make_model()
    stats = {}
    stats['orig_stmts'] = len(model.get_statements())
    stats['orig_monomers'] = len(pysb_model.monomers)
    stats['orig_rules'] = len(pysb_model.rules)
    # Extend the model with PMIDs
    stats['new_papers'], stats['new_abstracts'] = extend_model(model, pmids)
    stats['new_stmts'] = len(model.get_statements())
    pysb_model = model.make_model()
    stats['new_monomers'] = len(pysb_model.monomers)
    stats['new_rules'] = len(pysb_model.rules)
    # Save model
    model.save(inc_model_file)
    # Print and tweet the status message
    msg_str = make_status_message(stats)
    if msg_str is not None:
        print msg_str
        if use_twitter:
            twitter_client.update_status(msg_str, twitter_cred)
