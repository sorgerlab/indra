import re
import os
import sys
import imaplib
import email
import datetime
import getpass
import base64
import shutil

def get_mailboxes(M):
    res, mailboxes = M.list()
    if res == 'OK':
        return mailboxes
    else:
        return None

def select_mailbox(M, mailbox):
    res, data = M.select(mailbox)
    if res == 'OK':
        return data
    else:
        return None

def fetch_email(M, msg_id):
    res, data = M.fetch(msg_id, '(RFC822)')
    if res == 'OK':
        # Data here is a list with 1 element containing a tuple
        # whose 2nd element is a long string containing the email
        raw_msg_txt = data[0][1]
        msg = email.message_from_string(raw_msg_txt)
        return msg
    else:
        return None

def get_headers(msg):
    headers = {}
    for k in msg.keys():
        decode = email.Header.decode_header(msg[k])[0]
        txt = unicode(decode[0])
        headers[k] = txt
    return headers

def get_text(msg):
    parts = msg.get_payload()
    content_type = msg.get_content_type()
    msg_txt = None
    if content_type == 'text/html':
        if isinstance(parts, basestring):
            if parts.startswith('<?xml'):
                msg_txt = parts
            else:
                msg_txt = base64.b64decode(parts)
    else:
        print 'Can\'t handle content type %s' % content_type
    return msg_txt

def print_msg(msg):
    headers = get_headers(msg)
    text = get_text(msg)
    print '-----------'
    print 'Subject: %s' % headers['Subject']
    print 'From: %s' % headers['From']
    print 'To: %s' % headers['To']

    print 'Message:'
    print text

def get_message_pmids(M, day_limit=10):
    if day_limit is not None:
        date_now = datetime.datetime.now()
        date_rel = date_now - datetime.timedelta(days=10)
        date_str = date_rel.strftime('%d-%b-%Y')
        res, data = M.search(None, '(SINCE "%s")' % date_str)
    else:
        res, data = M.search(None, 'ALL')
    # Data here is a space-separated list of message IDs
    # like ['1 2 3']
    msg_ids = data[0].split(' ')
    pmids = []
    for mid in msg_ids:
        msg = fetch_email(M, mid)
        headers = get_headers(msg)
        subject = headers['Subject']
        subject_pmids = pmids_from_subject(subject)
        pmids += subject_pmids
        if headers['From'] == 'Sent by NCBI <nobody@ncbi.nlm.nih.gov>' or\
            headers['From'] == 'My NCBI <efback@ncbi.nlm.nih.gov>':
            text = get_text(msg)
            ncbi_pmids = pmids_from_ncbi_email(text)
            pmids += ncbi_pmids
    return pmids

def pmids_from_ncbi_email(msg_text):
    res = re.findall('PMID: [^.;]+', msg_text.replace('\n',''))
    pmids = [r[6:].strip() for r in res]
    return pmids

def pmids_from_subject(subject):
    pmids = []
    # TODO: this only works if the subject has PMIDxxx as a word
    # separated by spaces from other text.
    # We should use regexp to isolate the PMID
    subject_words = subject.split(' ')
    for w in subject_words:
        if w.startswith('PMID'):
            pmids.append(w[4:])
    return pmids

def gmail_login(email_addr, passwd):
    M = imaplib.IMAP4_SSL('imap.gmail.com')
    try:
        M.login(email_addr, passwd)
    except imaplib.IMAP4.error:
        print 'Login failed'
        return None
    return M

