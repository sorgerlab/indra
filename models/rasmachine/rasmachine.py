import os
import time
import shutil
import getpass
import gmail_client
import twitter_client
from incremental_model import IncrementalModel
from indra import reach
from indra.literature import pubmed_client, get_full_text

model_path = os.path.dirname(os.path.abspath(__file__))
use_twitter = False

def get_email_pmids(uname, passfile):
    try:
        passwd = open(pass_file, 'rt').read().strip()
    except IOError:
        passwd = getpass.getpass()

    M = gmail_client.gmail_login(uname, passwd)
    gmail_client.select_mailbox(M, 'INBOX')

    pmids = gmail_client.get_message_pmids(M)
    print 'Collected %d PMIDs' % len(pmids)
    return pmids

def process_paper(pmid):
    abstract_path = model_path + '/jsons/abstract/PMID%s.json' % pmid
    fulltext_path = model_path + '/jsons/full/PMID%s.json' % pmid

    # If the paper has been parsed, use the parse output file
    if os.path.exists(abstract_path):
        rp = reach.process_json_file(abstract_path, citation=pmid)
        txt_format = 'abstract'
    elif os.path.exists(fulltext_path):
        rp = reach.process_json_file(fulltext_path, citation=pmid)
        txt_format = 'txt'
    # If the paper has not been parsed, download the text and parse
    else:
        txt, txt_format = indra.literature.get_full_text(pmid)
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
    print time.strftime('%c')

    # Get email PMIDs
    uname = 'therasmachine@gmail.com'
    pass_file = 'rasmachine_cred.txt'
    pmids = get_email_pmids(uname, pass_file)
    # Load the model
    rasmodel = IncrementalModel(model_path + '/rasmodel.pkl')
    pysb_model = rasmodel.make_model()
    stats = {}
    stats['orig_stmts'] = len(rasmodel.get_statements())
    stats['orig_monomers'] = len(pysb_model.monomers)
    stats['orig_rules'] = len(pysb_model.rules)
    # Extend the model with PMIDs
    stats['new_papers'], stats['new_abstracts'] = extend_model(rasmodel, pmids)
    stats['new_stmts'] = len(rasmodel.get_statements())
    pysb_model = rasmodel.make_model()
    stats['new_monomers'] = len(pysb_model.monomers)
    stats['new_rules'] = len(pysb_model.rules)
    # Save model
    rasmodel.save(model_path + '/rasmodel.pkl')
    # Print and tweet the status message
    msg_str = make_status_message(stats)
    if msg_str is not None:
        print msg_str
        if use_twitter:
            twitter_client.update_status(msg_str, 'twitter_cred.txt')
    print model
