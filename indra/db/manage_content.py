from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, int

from sys import version_info, exit
if version_info.major is not 3:
    msg = "Python 3.x is required to use this module."
    if __name__ ==  '__main__':
        print(msg)
        exit()
    else:
        raise ImportError(msg)

import csv
import time
import tarfile
import zlib
import logging
import xml.etree.ElementTree as ET
import re
import os
import multiprocessing as mp
import pickle
import sys
from os import path
from ftplib import FTP
from io import BytesIO

logger = logging.getLogger('content_manager')
if __name__ == '__main__':
    # NOTE: PEP8 will complain about this, however having the args parsed up
    # here prevents a long wait just to fined out you entered a command wrong.
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description='Manage content on INDRA\'s database.'
        )
    parser.add_argument(
        choices=['upload', 'update'],
        dest='task',
        help=('Choose whether you want to perform an initial upload or update '
              'the existing content on the database.')
        )
    parser.add_argument(
        '-c', '--continue',
        dest='continuing',
        action='store_true',
        help='Continue uploading or updating, picking up where you left off.'
        )
    parser.add_argument(
        '-n', '--num_procs',
        dest='num_procs',
        type=int,
        default=1,
        help=('Select the number of processors to use during this operation. '
              'Default is 1.')
        )
    parser.add_argument(
        '-d', '--debug',
        dest='debug',
        action='store_true',
        help='Run with debugging level output.'
        )
    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
        from indra.db import logger as db_logger
        db_logger.setLevel(logging.DEBUG)

    if not args.continuing and args.task == 'upload':
        print("#"*63)
        print(
            "# You are about to wipe the database completely clean and     #\n"
            "# load it from scratch, which could take hours, and eliminate #\n"
            "# any updates, and potentially interfere with other users.    #\n"
            "#                                                             #\n"
            "# If you wish to continue an earlier upload, add the -c       #\n"
            "# option, and if you wish to update the database, specify     #\n"
            "# `update` in your command. For mor details, use --help.      #"
            )
        print("#"*63)
        resp = input("Are you sure you want to continue? [yes/NO]: ")
        if resp == 'no':
            print ("Aborting...")
            exit()
from indra.util import zip_string
from indra.util import UnicodeXMLTreeBuilder as UTB
from indra.literature.pmc_client import id_lookup
from indra.literature import pubmed_client
from indra.db import get_primary_db, texttypes, formats
from indra.db import sql_expressions as sql_exp


try:
    from psycopg2 import IntegrityError
except ImportError:
    class IntegrityError(object):
        "Using this in a try-except will catch nothing. (That's the point.)"
        pass


ftp_blocksize = 33554432  # Chunk size recommended by NCBI
BATCH_SIZE = 10000


class UploadError(Exception):
    pass


class NihFtpClient(object):
    """High level access to the NIH FTP repositories.

    Parameters
    ----------
    my_path : str
        The path to the subdirectory around which this client operates.
    ftp_url : str
        The url to the ftp site. May be a local directory (see `local`). By
        default this is `'ftp.ncbi.nlm.nih.gov'`.
    local : bool
        These methods may be run on a local directory (intended for testing).
        (default is `False`).
    """
    def __init__(self, my_path, ftp_url='ftp.ncbi.nlm.nih.gov', local=False):
        self.my_path = my_path
        self.is_local = local
        self.ftp_url = ftp_url

    def _path_join(self, *args):
        joined_str = path.join(*args)
        part_list = joined_str.split('/')
        for part in part_list[1:]:
            if part == '..':
                idx = part_list.index(part) - 1
                part_list.pop(idx)
                part_list.pop(idx)
        ret = path.join(*part_list)
        if part_list[0] == '':
            ret = '/' + ret
        return ret

    def get_ftp_connection(self, ftp_path=None):
        if ftp_path is None:
            ftp_path = self.my_path
        # Get an FTP connection
        ftp = FTP(self.ftp_url)
        ftp.login()
        # Change to the manuscripts directory
        ftp.cwd(ftp_path)
        return ftp

    def get_xml_file(self, xml_file):
        "Get the content from an xml file as an ElementTree."
        logger.info("Downloading %s" % (xml_file))
        xml_bytes = self.get_uncompressed_bytes(xml_file, force_str=False)
        logger.info("Parsing XML metadata")
        return ET.XML(xml_bytes, parser=UTB())

    def get_csv_as_dict(self, csv_file, cols=None, header=None):
        "Get the content from a csv file as a list of dicts."
        csv_str = self.get_file(csv_file)
        lst = [row for row in csv.reader(csv_str.splitlines())]
        if cols is None:
            if header is not None:
                cols = lst[header]
                lst = lst[header+1:]
            else:
                cols = list(range(len(lst[0])))
        return [dict(zip(cols, row)) for row in lst]

    def ret_file(self, f_path, buf):
        "Load the content of a file into the given buffer."
        full_path = self._path_join(self.my_path, f_path)
        if not self.is_local:
            with self.get_ftp_connection() as ftp:
                ftp.retrbinary('RETR /%s' % full_path,
                               callback=lambda s: buf.write(s),
                               blocksize=ftp_blocksize)
                buf.flush()
        else:
            with open(self._path_join(self.ftp_url, full_path), 'rb') as f:
                buf.write(f.read())
                buf.flush()
        return

    def download_file(self, f_path):
        "Download a file into a file given by f_path."
        name = path.basename(f_path)
        with open(name, 'wb') as gzf:
            self.ret_file(f_path, gzf)
        return name

    def get_file(self, f_path, force_str=True):
        "Get the contents of a file as a string."
        gzf_bytes = BytesIO()
        self.ret_file(f_path, gzf_bytes)
        ret = gzf_bytes.getvalue()
        if force_str and isinstance(ret, bytes):
            ret = ret.decode('utf8')
        return ret

    def get_uncompressed_bytes(self, f_path, force_str=True):
        "Get a file that is gzipped, and return the unzipped string."
        zipped_str = self.get_file(f_path, force_str=False)
        str_bytes = zlib.decompress(zipped_str, 16+zlib.MAX_WBITS)
        if force_str:
            ret = str_bytes.decode('utf8')
        else:
            ret = str_bytes
        return ret

    def ftp_ls(self, ftp_path=None):
        "Get a list of the contents in the ftp directory."
        if ftp_path is None:
            ftp_path = self.my_path
        else:
            ftp_path = self._path_join(self.my_path, ftp_path)
        if not self.is_local:
            with self.get_ftp_connection() as ftp:
                contents = ftp.nlst()
        else:
            contents = os.listdir(self._path_join(self.ftp_url, ftp_path))
        return contents


class Manager(object):
    """Abstract class for all upload/update managers.

    This abstract class provides the api required for any object that is
    used to manage content between the database and the content.
    """
    my_source = NotImplemented

    def copy_into_db(self, db, *args, **kwargs):
        "Wrapper around the db.copy feature, pickels args upon exception."
        try:
            db.copy(*args, **kwargs)
        except IntegrityError as e:
            pkl_file_fmt = "copy_failure_%d.pkl"
            i = 0
            while os.path.exists(pkl_file_fmt % i):
                i += 1
            with open(pkl_file_fmt % i, 'wb') as f:
                pickle.dump((e, args, kwargs), f, protocol=3)
            logger.error('Failed in a copy. Pickling args and kwargs.')
            logger.exception(e)
            logger.info('Continuing...')

    def populate(self, db):
        "A stub for the method used to initially populate the database."
        raise NotImplementedError(
            "`Populate` not implemented for `%s`." % self.__class__.__name__
            )

    def update(self, db):
        "A stub for the method used to update the content on the database."
        raise NotImplementedError(
            "`Update` not implemented for `%s`." % self.__class__.__name__
            )


class NihManager(Manager):
    """Abstract class for all the managers that use the NIH FTP service.

    See `NihFtpClient` for parameters.
    """
    my_path = NotImplemented

    def __init__(self, *args, **kwargs):
        self.ftp = NihFtpClient(self.my_path, *args, **kwargs)


class Medline(NihManager):
    "Manager for the medline content."
    my_path = 'pubmed/baseline'
    my_source = 'pubmed'

    def get_deleted_pmids(self):
        del_pmid_str = self.ftp.get_uncompressed_bytes(
            '../deleted.pmids.gz'
            )
        pmid_list = [
            line.strip() for line in del_pmid_str.split('\n')
            ]
        return pmid_list

    def get_file_list(self):
        all_files = self.ftp.ftp_ls()
        return [k for k in all_files if k.endswith('.xml.gz')]

    def get_article_info(self, xml_file, q=None):
        tree = self.ftp.get_xml_file(xml_file)
        article_info = pubmed_client.get_metadata_from_xml_tree(
            tree,
            get_abstracts=True,
            prepend_title=False
            )
        if q is not None:
            q.put((xml_file, article_info))
            return
        else:
            return article_info

    def fix_doi(self, doi):
        "Sometimes the doi is doubled (no idea why). Fix it."
        if doi is None:
            return doi
        L = len(doi)
        if L % 2 is not 0:
            return
        if doi[:L] != doi[L:]:
            return
        logger.info("Fixing doubled doi: %s" % doi)
        return doi[:L]

    def upload_article(self, db, article_info):
        deleted_pmids = self.get_deleted_pmids()

        logger.info("%d PMIDs in XML dataset" % len(article_info))
        # Convert the article_info into a list of tuples for insertion into
        # the text_ref table
        text_ref_records = []
        text_content_info = {}
        valid_pmids = set(article_info.keys()).difference(set(deleted_pmids))
        logger.info("%d valid PMIDs" % len(valid_pmids))
        existing_pmids = set(db.get_pmids(valid_pmids))
        logger.info(
            "%d valid PMIDs already in text_refs." % len(existing_pmids)
            )
        pmids_to_add = valid_pmids.difference(existing_pmids)
        logger.info("%d PMIDs to add to text_refs" % len(pmids_to_add))
        for pmid in pmids_to_add:
            pmid_data = article_info[pmid]
            rec = (
                pmid, pmid_data.get('pmcid'),
                self.fix_doi(pmid_data.get('doi')),
                pmid_data.get('pii')
                )
            text_ref_records.append(
                tuple([None if not r else r for r in rec])
                )
            abstract = pmid_data.get('abstract')
            # Make sure it's not an empty or whitespace-only string
            if abstract and abstract.strip():
                abstract_gz = zip_string(abstract)
                text_content_info[pmid] = (self.my_source, formats.TEXT,
                                           texttypes.ABSTRACT, abstract_gz)

        self.copy_into_db(
            db,
            'text_ref',
            text_ref_records,
            ('pmid', 'pmcid', 'doi', 'pii',)
            )

        # Build a dict mapping PMIDs to text_ref IDs
        pmid_list = list(text_content_info.keys())
        tref_list = db.select_all(
            'text_ref',
            db.TextRef.pmid.in_([p for p in pmid_list])
            )
        pmid_tr_dict = {pmid: trid for (pmid, trid) in
                        db.get_values(tref_list, ['pmid', 'id'])}

        # Add the text_ref IDs to the content to be inserted
        text_content_records = []
        for pmid, tc_data in text_content_info.items():
            if pmid not in pmid_tr_dict.keys():
                continue
            tr_id = pmid_tr_dict[pmid]
            text_content_records.append((tr_id,) + tc_data)

        self.copy_into_db(
            db,
            'text_content',
            text_content_records,
            cols=('text_ref_id', 'source', 'format', 'text_type', 'content',)
            )
        return True

    def populate(self, db, n_procs=1, continuing=False):
        """Perform the initial input of the pubmed content into the database.

        Parameters
        ----------
        db : indra.db.DatabaseManager instance
            The database to which the data will be uploaded.
        n_procs : int
            The number of processes to use when parsing xmls.
        continuing : bool
            If true, assume that we are picking up after an error, or otherwise
            continuing from an earlier process. This means we will skip over
            source files contained in the database. If false, all files will be
            read and parsed.
        """
        xml_files = self.get_file_list()
        sf_list = db.select_all(
            'source_file',
            db.SourceFile.source == self.my_source
            )
        existing_files = [sf.name for sf in sf_list]

        # This could perhaps be simplified with map_async from mp.pool.
        q = mp.Queue()
        proc_list = []
        for xml_file in xml_files:
            if continuing and xml_file in existing_files:
                logger.info("Skipping %s. Already uploaded." % xml_file)
                continue
            p = mp.Process(
                target=self.get_article_info,
                args=(xml_file, q, ),
                daemon=True
                )
            proc_list.append(p)
        n_tot = len(proc_list)

        for _ in range(n_procs):
            if len(proc_list):
                proc_list.pop(0).start()

        while len(proc_list):
            xml_file, article_info = q.get()  # Block until at least 1 is done.
            proc_list.pop(0).start()
            self.upload_article(db, article_info)
            if xml_file not in existing_files:
                db.insert('source_file', source=self.my_source, name=xml_file)
            n_tot -= 1

        while n_tot is not 0:
            xml_file, article_info = q.get()
            if xml_file not in existing_files:
                db.insert('source_file', source=self.my_source, name=xml_file)
            self.upload_article(db, article_info)
            n_tot -= 1

        return


class PmcManager(NihManager):
    """Abstract class for uploaders of PMC content.

    For Paramters, see `NihManager`.
    """
    my_source = NotImplemented

    def __init__(self, *args, **kwargs):
        super(PmcManager, self).__init__(*args, **kwargs)
        self.tr_cols = ('pmid', 'pmcid', 'doi', 'manuscript_id',)
        self.tc_cols = ('text_ref_id', 'source', 'format', 'text_type',
                        'content',)

    def get_missing_pmids(self, tr_data):
        "Try to get missing pmids using the pmc client."
        num_missing = 0
        num_found = 0

        logger.debug("Getting missing pmids.")

        # TODO: This is very slow...should find a way to speed it up.
        for tr_entry in tr_data:
            if tr_entry['pmid'] is None:
                num_missing += 1
                ret = id_lookup(tr_entry['pmcid'])
                if 'pmid' in ret.keys():
                    tr_entry['pmid'] = ret['pmid']
                    num_found += 1

        ''' # The web api does not support this much access, sadly.
        thread_list = []
        for tr_entry in tr_data:
            if tr_entry['pmid'] is None:
                th = Thread(target=lookup_pmid, args=[tr_entry])
                thread_list.append(th)

        N = min(10, len(thread_list))
        logger.debug("Starting %d threading pool." % N)
        active_threads = []
        for _ in range(N):
            th = thread_list.pop()
            th.start()
            active_threads.append(th)

        while len(thread_list):
            for th in active_threads[:]:
                if not th.is_alive():
                    th.join()
                    active_threads.remove(th)
                    if len(thread_list):
                        new_th = thread_list.pop()
                        new_th.start()
                        active_threads.append(th)
            sleep(0.1)

        for th in active_threads:
            th.join()
        '''
        logger.debug("Found %d/%d new pmids." % (num_found, num_missing))
        return

    def filter_text_refs(self, db, tr_data):
        "Try to reconcile the data we have with what's already on the db."
        logger.info("Beginning to filter text refs...")

        # This is a helper for accessing the data tuples we create
        def id_idx(id_type):
            return self.tr_cols.index(id_type)

        # If there are not actual refs to work with, don't waste time.
        if not len(tr_data):
            return [], []

        # Check for any pmids we can get from the pmc client (this is slow!)
        self.get_missing_pmids(tr_data)

        # Turn the list of dicts into a set of tuples
        tr_data_set = {tuple([entry[id_type] for id_type in self.tr_cols])
                       for entry in tr_data}

        # Get all text refs that match any of the id data we have. This means
        # each text ref WILL find a match in the data we have (unless something
        # is seriously broken.
        or_list = []
        for id_type in self.tr_cols:
            id_list = [entry[id_type] for entry in tr_data
                       if entry[id_type] is not None]
            if id_list:
                or_list.append(getattr(db.TextRef, id_type).in_(id_list))
        tr_list = db.select_all(db.TextRef, sql_exp.or_(*or_list))

        # Create an index of tupled data entries for quick lookups by any id
        # type, for example tr_data_idx_dict['pmid'][<a pmid>] will get the
        # tuple with all the id data. This avoids several iterations through
        # the list of text ref data dicts on for each text ref, at the cost of
        # using extra memory.
        tr_data_idx_dict = {id_type: {e[id_idx(id_type)]: e
                                      for e in tr_data_set
                                      if e[id_idx(id_type)] is not None}
                            for id_type in self.tr_cols}

        # Look for updates to the existing text refs
        tr_data_matched_set = set()
        pmcids_to_skip = set()
        for tr in tr_list:
            match_set = set()

            # Find the matches in the data. We continue looking after finding
            # one match in case there are any inconsistencies which need to be
            # considered.
            for id_type, tr_data_idx in tr_data_idx_dict.items():
                candidate = tr_data_idx.get(getattr(tr, id_type))
                if candidate is not None:
                    match_set.add(candidate)

            # As per the process of getting the tr_list, every tr MUST have a
            # match, or else something is broken.
            assert match_set, "No matches found, which is impossible."

            # Assuming we found exactly one matched data entry, we can now look
            # for new id data (for example manuscript ids) and update the tr
            # objects. These changes are transmuted to the db via the commit
            # command below.
            if len(match_set) == 1:
                tr_new = match_set.pop()

                # This is how we tell what doesn't need to be added to the db.
                tr_data_matched_set.add(tr_new)

                # Go through all the id_types
                for i, id_type in enumerate(self.tr_cols):
                    # Check if the text ref is missing that id.
                    if getattr(tr, id_type) is None:
                        # If so, and if our new data does have that id, update
                        # the text ref.
                        if tr_new[i] is not None:
                            setattr(tr, id_type, tr_new[i])
                    else:
                        # Check to see that all the ids agree. If not, report
                        # it in the review.txt file.
                        # NOTE: If this is ever done on AWS or through a
                        # container, the review file MUST be loaded somewhere
                        # it won't disappear. (such as s3). Perhaps these could
                        # be logged on the database?
                        if tr_new[i] is not None \
                         and tr_new[i] != getattr(tr, id_type):
                            with open('review.txt', 'a') as f:
                                f.write(
                                    'Got conflicting id data: in db %s vs %s.'
                                    % ([getattr(tr, id_type)
                                        for id_type in self.tr_cols],
                                       [tr_new[i]
                                        for i in range(len(self.tr_cols))])
                                    )
                            # If the conflict was with a pmcid, don't try to
                            # add the text content.
                            if id_type == 'pmcid':
                                pmcids_to_skip.add(
                                    tr_new[id_idx('pmcid')]
                                    )
            else:
                # These still matched something in the db, so they shouldn't be
                # uploaded as new refs.
                for tr_new in match_set:
                    tr_data_matched_set.add(tr_new)
                    pmcids_to_skip.add(tr_new[id_idx('pmcid')])

                # This condition only occurs if the records we got are
                # internally inconsistent. This is rare, but it can happen.
                with open('review.txt', 'a') as f:
                    f.write(
                        ('Got multiple matches for data from %s: %s.'
                         ' Please review.\n')
                        % (self.my_source, match_set)
                        )

        # This applies all the changes made to the text refs to the db.
        db.commit("Failed to update with new ids.")

        # Now update the text refs with any new refs that were found
        filtered_tr_records = tr_data_set.difference(tr_data_matched_set)

        return filtered_tr_records, pmcids_to_skip

    def filter_text_content(self, db, tc_data):
        'Filter the text content to identify pre-existing records.'
        logger.info("Beginning to filter text content...")
        arc_pmcid_list = [tc['pmcid'] for tc in tc_data]
        if not len(tc_data):
            return []
        tref_list = db.select_all(
            'text_ref',
            db.TextRef.pmcid.in_(arc_pmcid_list)
            )
        pmcid_trid_dict = {
            pmcid: trid for (pmcid, trid) in
            db.get_values(tref_list, ['pmcid', 'id'])
            }
        # This should be a very small list, in general.
        existing_tcs = db.select_all(
            'text_content',
            db.TextContent.text_ref_id.in_(pmcid_trid_dict.values()),
            db.TextContent.source == self.my_source,
            db.TextContent.format == formats.XML
            )
        existing_tc_records = [
            (tc.text_ref_id, tc.source, tc.format, tc.text_type)
            for tc in existing_tcs
            ]
        tc_records = []
        for tc in tc_data:
            tc_records.append(
                (
                    pmcid_trid_dict[tc['pmcid']],
                    self.my_source,
                    formats.XML,
                    tc['text_type'],
                    tc['content']
                    )
                )
        filtered_tc_records = [
            rec for rec in tc_records if rec[:-1] not in existing_tc_records
            ]
        logger.info("Finished filtering the text content...")
        return list(set(filtered_tc_records))

    def upload_batch(self, db, tr_data, tc_data):
        "Add a batch of text refs and text content to the database."
        filtered_tr_records, pmcids_to_skip = \
            self.filter_text_refs(db, tr_data)
        if len(pmcids_to_skip) is not 0:
            mod_tc_data = [
                tc for tc in tc_data if tc['pmcid'] not in pmcids_to_skip
                ]
        else:
            mod_tc_data = tc_data

        # Upload the text content data.
        logger.info('Adding %d new text refs...' % len(filtered_tr_records))
        self.copy_into_db(
            db,
            'text_ref',
            filtered_tr_records,
            self.tr_cols
            )

        # Process the text content data
        filtered_tc_records = self.filter_text_content(db, mod_tc_data)

        # Upload the text content data.
        logger.info('Adding %d more text content entries...' %
                    len(filtered_tc_records))
        self.copy_into_db(
            db,
            'text_content',
            filtered_tc_records,
            self.tc_cols
            )

    def get_data_from_xml_str(self, xml_str, filename):
        "Get the data out of the xml string."
        try:
            tree = ET.XML(xml_str.encode('utf8'))
        except ET.ParseError:
            logger.info("Could not parse %s. Skipping." % filename)
            return None
        id_data = {
            e.get('pub-id-type'): e.text for e in
            tree.findall('.//article-id')
            }
        if 'pmc' not in id_data.keys():
            logger.info("Did not get a 'pmc' in %s." % filename)
            return None
        if 'pmcid' not in id_data.keys():
            id_data['pmcid'] = 'PMC' + id_data['pmc']
        if 'manuscript' in id_data.keys():
            id_data['manuscript_id'] = id_data['manuscript']
        tr_datum = {k: id_data.get(k) for k in self.tr_cols}
        tc_datum = {
            'pmcid': id_data['pmcid'],
            'text_type': texttypes.FULLTEXT,
            'content': zip_string(xml_str)
            }
        return tr_datum, tc_datum

    def upload_archive(self, archive, q=None, db=None):
        "Process a single tar gzipped article archive."
        tr_data = []
        tc_data = []

        def submit(tag, tr_data, tc_data):
            batch_name = 'final batch' if tag is 'final' else 'batch %d' % tag
            logger.info("Submitting %s of data for %s..." %
                        (batch_name, archive))

            if q is not None:
                q.put(((batch_name, archive), tr_data[:], tc_data[:]))
            elif db is not None:
                self.upload_batch(db, tr_data[:], tc_data[:])
            else:
                raise UploadError(
                    "upload_archive must receive either a db instance"
                    " or a queue instance."
                    )
            tr_data.clear()
            tc_data.clear()

        with tarfile.open(archive, mode='r:gz') as tar:
            logger.info('Loading %s...' % archive)
            xml_files = [m for m in tar.getmembers() if m.isfile()]
            for i, xml_file in enumerate(xml_files):
                xml_str = tar.extractfile(xml_file).read().decode('utf8')
                res = self.get_data_from_xml_str(xml_str, xml_file.name)
                if res is None:
                    continue
                else:
                    tr, tc = res
                tr_data.append(tr)
                tc_data.append(tc)
                if (i+1) % BATCH_SIZE is 0:
                    submit((i+1)/BATCH_SIZE, tr_data, tc_data)
            else:
                submit('final', tr_data, tc_data)
        return

    def process_archive(self, archive, q=None, db=None):
        try:
            logger.info('Downloading archive %s.' % archive)
            self.ftp.download_file(archive)
            self.upload_archive(archive, q=q, db=db)
        finally:
            os.remove(archive)

    def get_file_list(self):
        return [k for k in self.ftp.ftp_ls() if self.is_archive(k)]

    def populate(self, db, n_procs=1, continuing=False):
        """Perform the initial population of the pmc content into the database.

        Parameters
        ----------
        db : indra.db.DatabaseManager instance
            The database to which the data will be uploaded.
        n_procs : int
            The number of processes to use when parsing xmls.
        continuing : bool
            If true, assume that we are picking up after an error, or
            otherwise continuing from an earlier process. This means we will
            skip over source files contained in the database. If false, all
            files will be read and parsed.
        """
        archives = self.get_file_list()

        sf_list = db.select_all(
            'source_file',
            db.SourceFile.source == self.my_source
            )
        existing_arcs = [sf.name for sf in sf_list]

        q = mp.Queue(len(archives))
        wait_list = []
        for archive in archives:
            if continuing and archive in existing_arcs:
                logger.info("Skipping %s. Already uploaded." % archive)
                continue
            p = mp.Process(
                target=self.process_archive,
                args=(archive, ),
                kwargs={'q': q, },
                daemon=True
                )
            wait_list.append((archive, p))

        active_list = []

        def start_next_proc():
            if len(wait_list) is not 0:
                archive, proc = wait_list.pop(0)
                proc.start()
                active_list.append((archive, proc))

        # Start the processes running
        for _ in range(min(n_procs, len(archives))):
            start_next_proc()

        # Monitor the processes while any are still active.
        while len(active_list) is not 0:
            for a, p in [(a, p) for a, p in active_list if not p.is_alive()]:
                if a not in existing_arcs:
                    db.insert('source_file', source=self.my_source, name=a)
                active_list.remove((a, p))
                start_next_proc()
            try:
                # This will not block until at least one is done
                label, tr_data, tc_data = q.get_nowait()
            except Exception:
                continue
            logger.info("Beginning to upload %s from %s..." % label)
            self.upload_batch(db, tr_data, tc_data)
            logger.info("Finished %s from %s..." % label)
            time.sleep(0.1)

        # Empty the queue.
        while not q.empty():
            try:
                tr_data, tc_data = q.get(timeout=1)
            except Exception:
                break
            self.upload_batch(db, tr_data, tc_data)

        return


class PmcOA(PmcManager):
    "Manager for the pmc open access content."
    my_path = 'pub/pmc'
    my_source = 'pmc_oa'

    def is_archive(self, k):
        return k.startswith('articles') and k.endswith('.xml.tar.gz')


class Manuscripts(PmcManager):
    "Manager for the pmc manuscripts."
    my_path = 'pub/pmc/manuscript'
    my_source = 'manuscripts'

    def get_tarname_from_filename(self, fname):
        "Get the name of the tar file based on the file name (or a pmcid)."
        re_match = re.match('(PMC00\d).*?', fname)
        if re_match is not None:
            tarname = re_match.group(0) + 6*'X' + '.xml.tar.gz'
        else:
            tarname = None
        return tarname

    def is_archive(self, k):
        return k.endswith('.xml.tar.gz')

    def enrich_textrefs(self, db):
        """Method to add manuscript_ids to the text refs."""
        tr_list = db.select_all(db.TextRef,
                                db.TextContent.text_ref_id == db.TextRef.id,
                                db.TextContent.source == self.my_source,
                                db.TextRef.manuscript_id.is_(None))
        file_list = self.ftp.get_csv_as_dict('filelist.csv', header=0)
        pmcid_mid_dict = {entry['PMCID']: entry['MID'] for entry in file_list}
        pmid_mid_dict = {entry['PMID']: entry['MID'] for entry in file_list
                         if entry['PMID'] is not '0'}
        for tr in tr_list:
            if tr.pmcid is not None:
                tr.manuscript_id = pmcid_mid_dict[tr.pmcid]
            elif tr.pmid is not None and tr.pmid in pmid_mid_dict.keys():
                tr.manuscript_id = pmid_mid_dict[tr.pmid]
        db.commit("Could not update text refs with manuscript ids.")
        return


if __name__ == '__main__':
    db = get_primary_db()
    logger.info("Performing %s." % args.task)
    if args.task == 'upload':
        if not args.continuing:
            logger.info("Clearing TextContent and TextRef tables.")
            clear_succeeded = db._clear([db.TextContent, db.TextRef, db.SourceFile])
            if not clear_succeeded:
                sys.exit()
        Medline().populate(db, args.num_procs, args.continuing)
        PmcOA().populate(db, args.num_procs, args.continuing)
        Manuscripts().populate(db, args.num_procs, args.continuing)
    elif args.task == 'update':
        logger.warning("Sorry, this feature not yet available.")

    # High-level content update procedure
    # 1. Download MEDLINE baseline, will contain all PMIDs, abstracts,
    #    other info
    # 2. Download PMC, update text refs with PMCIDs where possible, with
    #    update commands.
    # 3. Add PMC-OA content.
    # 4. Add PMC-author manuscript information, updating files with manuscript
    #    IDs and content.
    # 5. (Optional): Run script to check for/obtain DOIs for those that are
    #    missing
    # 6. Obtain Elsevier XML for Elsevier articles
    # 7. Obtain Springer content for Springer articles
    # 8. Obtain scraped PDF content for other articles available via OA-DOI.
