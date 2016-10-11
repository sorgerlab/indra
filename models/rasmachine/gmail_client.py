from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, bytes
import re
import os
import sys
import imaplib
import email
import email.header
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
    """Returns the given email message as a unicode string."""
    res, data = M.fetch(msg_id, '(RFC822)')
    if res == 'OK':
        # Data here is a list with 1 element containing a tuple
        # whose 2nd element is a long string containing the email
        # The content is a bytes that must be decoded
        raw_msg_txt = data[0][1]
        # In Python3, we call message_from_bytes, but this function doesn't
        # exist in Python 2.
        try:
            msg = email.message_from_bytes(raw_msg_txt)
        except AttributeError:
            msg = email.message_from_string(raw_msg_txt)
        # At this point, we have a message containing bytes (not unicode)
        # fields that will still need to be decoded, ideally according to the
        # character set specified in the message.
        return msg
    else:
        return None

def get_headers(msg):
    """Takes email.message.Message object initialized from unicode string,
    returns dict with header fields."""
    headers = {}
    for k in msg.keys():
        # decode_header decodes header but does not convert charset, so these
        # may still be bytes, even in Python 3. However, if it's ASCII
        # only (hence unambiguous encoding), the header fields come back
        # as str (unicode) in Python 3.
        (header_txt, charset) = email.header.decode_header(msg[k])[0]
        if charset is not None:
            header_txt = header_txt.decode(charset)
        headers[k] = header_txt
    return headers

def get_text(msg):
    # Decode=True argument handles quoted-printable and Base64 encoding.
    # parts variable may be a string or a list (in the case of a multi-part
    # message).
    parts = msg.get_payload(decode=True)
    content_type = msg.get_content_type()
    msg_txt = None
    if content_type == 'text/html':
        if isinstance(parts, list) or isinstance(parts, tuple):
            pass
        # If this is a bytes, we need to decode it
        elif isinstance(parts, bytes):
            charset = msg.get_charset()
            if charset is None:
                msg_txt = parts.decode('utf-8')
            else:
                msg_txt = parts.decode(charset)
        # If it's already a str, we're good to go
        elif isinstance(parts, str):
            msg_txt = parts
        else:
            raise Exception("Message payload was neither string nor list.")
    else:
        print('Can\'t handle content type %s' % content_type)
    return msg_txt

def print_msg(msg):
    headers = get_headers(msg)
    text = get_text(msg)
    print('-----------')
    print('Subject: %s' % headers['Subject'])
    print('From: %s' % headers['From'])
    print('To: %s' % headers['To'])

    print('Message:')
    print(text)

def get_message_pmids(M, day_limit=10):
    if day_limit is not None:
        date_now = datetime.datetime.now()
        date_rel = date_now - datetime.timedelta(days=10) # 10
        date_str = date_rel.strftime('%d-%b-%Y')
        res, data = M.search(None, '(SINCE "%s")' % date_str)
    else:
        res, data = M.search(None, 'ALL')
    # Data here is a space-separated list of message IDs
    # like ['1 2 3']
    msg_ids_str = data[0].decode('utf-8')
    msg_ids = msg_ids_str.split(' ')
    pmids = []
    for mid in msg_ids:
        # msg is returned as object containing bytes
        msg = fetch_email(M, mid)
        # get_headers converts fields to unicode
        headers = get_headers(msg)
        subject = headers['Subject']
        subject_pmids = pmids_from_subject(subject)
        pmids += subject_pmids
        if headers['From'] == 'Sent by NCBI <nobody@ncbi.nlm.nih.gov>' or\
            headers['From'] == 'My NCBI <efback@ncbi.nlm.nih.gov>':
            # Returns unicode
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
        print('Login failed')
        return None
    return M

