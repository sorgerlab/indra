import json
import boto3
import logging
from os import environ
from indra.statements import stmts_to_json, stmts_from_json
from indra.statements.io import stmts_to_json_file, stmts_from_json_file
from . import file_defaults, InvalidCorpusError, CACHE, default_bucket, \
    default_key_base, default_profile
from .util import _stmts_dict_to_json, _json_to_stmts_dict, _json_dumper, \
    _json_loader

logger = logging.getLogger(__name__)


class Corpus(object):
    """Represent a corpus of statements with curation.

    Parameters
    ----------
    statements : list[indra.statement.Statement]
        A list of INDRA Statements to embed in the corpus.
    raw_statements : list[indra.statement.Statement]
        A List of raw statements forming the basis of the statements in
        'statements'.
    aws_name : str
        The name of the profile in the AWS credential file to use. 'default' is
        used by default.

    Attributes
    ----------
    statements : dict
        A dict of INDRA Statements keyed by UUID.
    raw_statements : list
        A list of the raw statements
    curations : dict
        A dict keeping track of the curations submitted so far for Statement
        UUIDs in the corpus.
    meta_data : dict
        A dict with meta data associated with the corpus
    """
    def __init__(self, corpus_id, statements=None, raw_statements=None,
                 meta_data=None, aws_name=default_profile):
        self.corpus_id = corpus_id
        self.statements = {st.uuid: st for st in statements} if statements \
            else {}
        self.raw_statements = raw_statements if raw_statements else []
        self.curations = {}
        self.meta_data = meta_data if meta_data else {}
        self.aws_name = aws_name
        self._s3 = None

    def _get_s3_client(self):
        if self._s3 is None:
            if environ.get('AWS_ACCESS_KEY_ID') and \
                    environ.get('AWS_SECRET_ACCESS_KEY'):
                logger.info('Got credentials in environment for client')
                self._s3 = boto3.session.Session(
                    aws_access_key_id=environ.get('AWS_ACCESS_KEY_ID'),
                    aws_secret_access_key=environ.get('AWS_SECRET_ACCESS_KEY')
                ).client('s3')
            else:
                logger.info('Using stored AWS profile for client')
                self._s3 = boto3.session.Session(
                    profile_name=self.aws_name).client('s3')
        return self._s3

    def __str__(self):
        return 'Corpus(%s -> %s)' % (str(self.statements), str(self.curations))

    def __repr__(self):
        return str(self)

    @classmethod
    def load_from_s3(cls, corpus_id, aws_name=default_profile,
                     bucket=default_bucket, force_s3_reload=False,
                     raise_exc=False):
        corpus = cls(corpus_id, statements=[], aws_name=aws_name)
        corpus.s3_get(bucket, cache=(not force_s3_reload),
                      raise_exc=raise_exc)
        return corpus

    def s3_put(self, bucket=default_bucket, cache=True):
        """Push a corpus object to S3 in the form of three json files

        The json files representing the object have S3 keys of the format
        <key_base_name>/<name>/<file>.json

        Parameters
        ----------
        bucket : str
            The S3 bucket to upload the Corpus to. Default: 'world-modelers'.
        cache : bool
            If True, also create a local cache of the corpus. Default: True.

        Returns
        -------
        keys : tuple(str)
            A tuple of three strings giving the S3 key to the pushed objects
        """
        # Note that the S3 path to each json file is of the form
        # <bucket>/indra_models/<corpus_id>/<file>.json"
        s3key = '%s/%s/' % (default_key_base, self.corpus_id)
        raw = s3key + file_defaults['raw'] + '.json'
        sts = s3key + file_defaults['sts'] + '.json'
        cur = s3key + file_defaults['cur'] + '.json'
        meta = s3key + file_defaults['meta'] + '.json'
        try:
            s3 = self._get_s3_client()
            # Structure and upload raw statements
            self._s3_put_file(s3,
                              raw,
                              json.dumps(stmts_to_json(self.raw_statements),
                                         indent=1),
                              bucket)

            # Structure and upload assembled statements
            self._s3_put_file(s3,
                              sts,
                              '\n'.join(json.dumps(jo, indent=1) for jo in
                                        _stmts_dict_to_json(self.statements)),
                              bucket)

            # Structure and upload curations
            self._s3_put_file(s3, cur, json.dumps(self.curations), bucket)

            # Upload meta data
            self._s3_put_file(s3, meta, json.dumps(self.meta_data), bucket)

            if cache:
                self._save_to_cache(raw, sts, cur)
            return list((raw, sts, cur))
        except Exception as e:
            logger.exception('Failed to put on s3: %s' % e)
            return None

    @staticmethod
    def _s3_put_file(s3, key, json_str, bucket=default_bucket):
        """Does the json.dumps operation for the the upload, i.e. json_obj
        must be an object that can be turned into a bytestring using
        json.dumps"""
        logger.info('Uploading %s to S3' % key)
        s3.put_object(Body=json_str, Bucket=bucket, Key=key)

    def _save_to_cache(self, raw=None, sts=None, cur=None, meta=None):
        """Helper method that saves the current state of the provided
        file keys"""
        # Assuming file keys are full s3 keys:
        # <base_name>/<dirname>/<file>.json

        # Raw:
        if raw:
            rawf = CACHE.joinpath(raw.replace(default_key_base + '/', ''))
            if not rawf.is_file():
                rawf.parent.mkdir(exist_ok=True, parents=True)
                rawf.touch(exist_ok=True)
            _json_dumper(jsonobj=stmts_to_json(self.raw_statements),
                         fpath=rawf.as_posix())

        # Assembled
        if sts:
            stsf = CACHE.joinpath(sts.replace(default_key_base + '/', ''))
            if not stsf.is_file():
                stsf.parent.mkdir(exist_ok=True, parents=True)
                stsf.touch(exist_ok=True)
            stmts_to_json_file(stmts=[s for _, s in self.statements.items()],
                               fname=stsf.as_posix(), format='jsonl')

        # Curation
        if cur:
            curf = CACHE.joinpath(cur.replace(default_key_base + '/', ''))
            if not curf.is_file():
                curf.parent.mkdir(exist_ok=True, parents=True)
                curf.touch(exist_ok=True)
            _json_dumper(jsonobj=self.curations, fpath=curf.as_posix())

        # Meta data
        if meta:
            metaf = CACHE.joinpath(meta.replace(default_key_base + '/', ''))
            if not metaf.is_file():
                metaf.parent.mkdir(exist_ok=True, parents=True)
                metaf.touch(exist_ok=True)
            _json_dumper(jsonobj=self.meta_data, fpath=metaf.as_posix())

    def s3_get(self, bucket=default_bucket, cache=True,
               raise_exc=False):
        """Fetch a corpus object from S3 in the form of three json files

        The json files representing the object have S3 keys of the format
        <s3key>/statements.json and <s3key>/raw_statements.json.

        Parameters
        ----------
        bucket : str
            The S3 bucket to fetch the Corpus from. Default: 'world-modelers'.
        cache : bool
            If True, look for corpus in local cache instead of loading it
            from s3. Default: True.
        raise_exc : bool
            If True, raise InvalidCorpusError when corpus failed to load

        """
        # Note that the S3 path to each json file is of the form
        # <bucket>/indra_models/<corpus_id>/<file>.json"
        s3key = '%s/%s/' % (default_key_base, self.corpus_id)
        raw = s3key + file_defaults['raw'] + '.json'
        sts = s3key + file_defaults['sts'] + '.json'
        cur = s3key + file_defaults['cur'] + '.json'
        meta = s3key + file_defaults['meta'] + '.json'
        try:
            logger.info('Loading corpus: %s' % s3key)
            s3 = self._get_s3_client()

            # Get and process raw statements
            raw_stmt_jsons = []
            if cache:
                raw_stmt_jsons = self._load_from_cache(raw) or []
            if not raw_stmt_jsons:
                raw_stmt_jsons_str = s3.get_object(
                    Bucket=bucket, Key=raw)['Body'].read()
                raw_stmt_jsons = json.loads(raw_stmt_jsons_str) or []
            self.raw_statements = stmts_from_json(raw_stmt_jsons)

            # Get and process assembled statements from list to dict
            json_stmts = []
            if cache:
                json_stmts = self._load_from_cache(sts) or []
            if not json_stmts:
                raw_str = s3.get_object(Bucket=bucket, Key=sts)[
                    'Body'].read().decode()
                if len(raw_str.split('\n')) > 1:
                    json_stmts = [json.loads(s) for s in raw_str.split('\n')]
                else:
                    json_stmts = json.loads(raw_str) or []

            self.statements = _json_to_stmts_dict(json_stmts)

            # Get and process curations if any
            curation_json = {}
            if cache:
                curation_json = self._load_from_cache(cur) or {}
            if not curation_json:
                curation_json = json.loads(s3.get_object(
                    Bucket=bucket, Key=cur)['Body'].read()) or {}
            self.curations = curation_json

            meta_json = {}
            try:
                if cache:
                    meta_json = self._load_from_cache(meta)
                if not meta_json:
                    meta_json = json.loads(s3.get_object(
                        Bucket=bucket, Key=meta)['Body'].read())
            except Exception as e:
                if isinstance(e, s3.exceptions.NoSuchKey):
                    logger.warning('No meta data found on s3')
                else:
                    logger.warning('No meta data found')
                meta_json = {}
            self.meta_data = meta_json

        except Exception as e:
            if raise_exc:
                raise InvalidCorpusError('Failed to get from s3: %s' % e)
            else:
                logger.warning('Failed to get from s3: %s' % e)

    def upload_curations(self, look_in_cache=False,
                         save_to_cache=False, bucket=default_bucket):
        """Upload the current state of curations for the corpus

        Parameters
        ----------
        look_in_cache : bool
            If True, when no curations are avaialbe check if there are
            curations cached locally. Default: False
        save_to_cache : bool
            If True, also save current curation state to cache. If
            look_in_cache is True, this option will have no effect. Default:
            False.
        bucket : str
            The bucket to upload to. Default: 'world-modelers'.
        """
        # Get curation file key
        file_key = '%s/%s/%s.json' % (default_key_base, self.corpus_id,
                                      file_defaults['cur'])
        # First see if we have any curations, then check in cache if
        # look_in_cache == True
        if self.curations:
            curations = self.curations
        elif look_in_cache:
            curations = self._load_from_cache(file_key)
        else:
            curations = None

        # Only upload if we actually have any curations to upload
        if curations:
            self._s3_put_file(s3=self._get_s3_client(),
                              key=file_key,
                              json_str=json.dumps(curations),
                              bucket=bucket)

        if self.curations and save_to_cache and not look_in_cache:
            self._save_to_cache(cur=file_key)

    def get_curations(self, look_in_cache=False):
        """Get curations for the corpus

        Parameters
        ----------
        look_in_cache : bool
            If True, look in local cache if there are no curations loaded

        Returns
        -------
        dict
            The curations for this corpus, if any
        """
        if self.curations:
            curations = self.curations
        elif look_in_cache:
            file_key = '%s/%s/%s.json' % (default_key_base, self.corpus_id,
                                          file_defaults['cur'])
            curations = self._load_from_cache(file_key) or {}
        else:
            curations = {}
        return curations

    @staticmethod
    def _load_from_cache(file_key):
        # Assuming file_key is cleaned, contains the file name and contains
        # the initial file base name:
        # <base_name>/<dirname>/<file>.json

        # Remove <base_name> and get local path to file
        local_file = CACHE.joinpath(
            '/'.join([s for s in file_key.split('/')[1:]]))

        # Load json object
        if local_file.is_file():
            if local_file.as_posix().endswith(file_defaults['sts'] + '.json'):
                return stmts_from_json_file(local_file.as_posix(),
                                            format='jsonl')
            else:
                return _json_loader(local_file.as_posix())
        return None

    def to_json_file(self, fname, w_newlines=False):
        """Dump the statements to a file in json format

        Parameters
        ----------
        fname : str
            A valid file path
        w_newlines : bool
            If True, the statements will be separated by newlines in the
            file. Default: False.
        """
        stmts_to_json_file(stmts=[s for _, s in self.statements.items()],
                           fname=fname,
                           format='jsonl' if w_newlines else 'json')
