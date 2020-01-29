"""This REST service allows real-time curation and belief updates for
a corpus of INDRA Statements."""
import json
import yaml
import boto3
import pickle
import logging
import argparse
from os import path
from pathlib import Path
from flask import Flask, request, jsonify, abort, Response
# Note: preserve EidosReader import as first one from indra
from indra.sources.eidos.reader import EidosReader
from indra.belief import BeliefEngine
from indra.tools import assemble_corpus as ac
from indra.belief.wm_scorer import get_eidos_bayesian_scorer
from indra.statements import stmts_from_json_file, stmts_to_json, \
    stmts_from_json, Statement
from indra.preassembler.hierarchy_manager import YamlHierarchyManager
from indra.preassembler.make_wm_ontologies import wm_ont_url, \
    load_yaml_from_url, rdf_graph_from_yaml


logger = logging.getLogger('live_curation')
app = Flask(__name__)
corpora = {}


default_bucket = 'world-modelers'
default_key_base = 'indra_models'
default_profile = 'wm'
file_defaults = {'raw': 'raw_statements',
                 'sts': 'statements',
                 'cur': 'curations'}

HERE = Path(path.abspath(__file__)).parent
CACHE = HERE.joinpath('_local_cache')
CACHE.mkdir(exist_ok=True)


def _json_loader(fpath):
    logger.info('Loading json file %s' % fpath)
    with open(fpath, 'r') as f:
        return json.load(f)


def _json_dumper(jsonobj, fpath):
    try:
        logger.info('Saving json object to file %s' % fpath)
        with open(fpath, 'w') as f:
            json.dump(obj=jsonobj, fp=f)
        return True
    except Exception as e:
        logger.error('Could not save json')
        logger.exception(e)
        return False


class Corpus(object):
    """Represent a corpus of statements with curation.

    Parameters
    ----------
    statements : list[indra.statement.Statement]
        A list of INDRA Statements to embed in the corpus.
    aws_name : str
        The name of the profile in the AWS credential file to use. 'default' is
        used by default.

    Attributes
    ----------
    statements : dict
        A dict of INDRA Statements keyed by UUID.
    curations : dict
        A dict keeping track of the curations submitted so far for Statement
        UUIDs in the corpus.
    """
    def __init__(self, statements, raw_statements=None,
                 aws_name=default_profile):
        self.statements = {st.uuid: st for st in statements}
        self.raw_statements = [] if not raw_statements else raw_statements
        self.curations = {}
        self.aws_name = aws_name
        self._s3 = None

    def _get_s3_client(self):
        if self._s3 is None:
            self._s3 = boto3.session.Session(
                profile_name=self.aws_name).client('s3')
        return self._s3

    def __str__(self):
        return 'Corpus(%s -> %s)' % (str(self.statements), str(self.curations))

    def __repr__(self):
        return str(self)

    @classmethod
    def load_from_s3(cls, s3key, aws_name=default_profile,
                     bucket=default_bucket, force_s3_reload=False,
                     raise_exc=False):
        corpus = cls([], aws_name=aws_name)
        corpus.s3_get(s3key, bucket, cache=(not force_s3_reload),
                      raise_exc=raise_exc)
        return corpus

    def s3_put(self, s3key, bucket=default_bucket, cache=True):
        """Push a corpus object to S3 in the form of three json files

        The json files representing the object have S3 keys of the format
        <key_base_name>/<name>/<file>.json

        Parameters
        ----------
        s3key : str
            The key to fetch the json files from. The key is assumed to be
            of the following form: "indra_models/<dirname>/<file>.json" and
            only <dirname> *must* be provided. Any combination of
            including/excluding 'indra_models' and/or <file>.json is
            accepted assuming the file ending '.json' is specified when
            <file>.json is specified.
        bucket : str
            The S3 bucket to upload the Corpus to. Default: 'world-modelers'.
        cache : bool
            If True, also create a local cache of the corpus. Default: True.

        Returns
        -------
        keys : tuple(str)
            A tuple of three strings giving the S3 key to the pushed objects
        """
        s3key = _clean_key(s3key) + '/'
        raw = s3key + file_defaults['raw'] + '.json'
        sts = s3key + file_defaults['sts'] + '.json'
        cur = s3key + file_defaults['cur'] + '.json'
        try:
            s3 = self._get_s3_client()
            # Structure and upload raw statements
            self._s3_put_file(s3, raw, stmts_to_json(self.raw_statements),
                              bucket)

            # Structure and upload assembled statements
            self._s3_put_file(s3, sts, _stmts_dict_to_json(self.statements),
                              bucket)

            # Structure and upload curations
            self._s3_put_file(s3, cur, self.curations, bucket)

            if cache:
                self._save_to_cache(raw, sts, cur)
            return list((raw, sts, cur))
        except Exception as e:
            logger.exception('Failed to put on s3: %s' % e)
            return None

    @staticmethod
    def _s3_put_file(s3, key, json_obj, bucket=default_bucket):
        """Does the json.dumps operation for the the upload, i.e. json_obj
        must be an object that can be turned into a bytestring using
        json.dumps"""
        logger.info('Uploading %s to S3' % key)
        s3.put_object(Body=json.dumps(json_obj),
                      Bucket=bucket, Key=key)

    def _save_to_cache(self, raw=None, sts=None, cur=None):
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
            _json_dumper(jsonobj=_stmts_dict_to_json(self.statements),
                         fpath=stsf.as_posix())

        # Curation
        if cur:
            curf = CACHE.joinpath(cur.replace(default_key_base + '/', ''))
            if not curf.is_file():
                curf.parent.mkdir(exist_ok=True, parents=True)
                curf.touch(exist_ok=True)
            _json_dumper(jsonobj=self.curations, fpath=curf.as_posix())

    def s3_get(self, s3key, bucket=default_bucket, cache=True,
               raise_exc=False):
        """Fetch a corpus object from S3 in the form of three json files

        The json files representing the object have S3 keys of the format
        <s3key>/statements.json and <s3key>/raw_statements.json.

        Parameters
        ----------
        s3key : str
            The key to fetch the json files from. The key is assumed to be
            of the following form: "indra_models/<dirname>/<file>.json" and
            only <dirname> *must* be provided. Any combination of
            including/excluding 'indra_models' and/or <file>.json is
            accepted assuming the file ending '.json' is specified when
            <file>.json is specified.
        bucket : str
            The S3 bucket to fetch the Corpus from. Default: 'world-modelers'.
        cache : bool
            If True, look for corpus in local cache instead of loading it
            from s3. Default: True.
        raise_exc : bool
            If True, raise InvalidCorpusError when corpus failed to load

        """
        s3key = _clean_key(s3key) + '/'
        raw = s3key + file_defaults['raw'] + '.json'
        sts = s3key + file_defaults['sts'] + '.json'
        cur = s3key + file_defaults['cur'] + '.json'
        try:
            logger.info('Loading corpus: %s' % s3key)
            s3 = self._get_s3_client()

            # Get and process raw statements
            raw_stmt_jsons = None
            if cache:
                raw_stmt_jsons = self._load_from_cache(raw)
            if raw_stmt_jsons is None:
                raw_stmt_jsons_str = s3.get_object(
                    Bucket=bucket, Key=raw)['Body'].read()
                raw_stmt_jsons = json.loads(raw_stmt_jsons_str)
            self.raw_statements = stmts_from_json(raw_stmt_jsons)

            # Get and process assembled statements from list to dict
            json_stmts = None
            if cache:
                json_stmts = self._load_from_cache(sts)
            if json_stmts is None:
                json_stmts = json.loads(s3.get_object(
                    Bucket=bucket, Key=sts)['Body'].read())

            self.statements = _json_to_stmts_dict(json_stmts)

            # Get and process curations if any
            curation_jsons = None
            if cache:
                curation_jsons = self._load_from_cache(cur)
            if curation_jsons is None:
                curation_jsons = json.loads(s3.get_object(
                    Bucket=bucket, Key=cur)['Body'].read())
            self.curations = {uid: c for uid, c in curation_jsons.items()}

        except Exception as e:
            if raise_exc:
                raise InvalidCorpusError('Failed to get from s3: %s' % e)
            else:
                logger.warning('Failed to get from s3: %s' % e)

    def upload_curations(self, corpus_id, look_in_cache=False,
                         save_to_cache=False, bucket=default_bucket):
        """Upload the current state of curations for the corpus

        Parameters
        ----------
        corpus_id : str
            The corpus ID of the curations to upload
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
        file_key = _clean_key(corpus_id) + '/' + \
                   file_defaults['cur'] + '.json'

        curations = self.curations if self.curations else (
            self._load_from_cache(file_key) if look_in_cache else None)

        self._s3_put_file(s3=self._get_s3_client(),
                          key=file_key,
                          json_obj=curations,
                          bucket=bucket)

        if save_to_cache and not look_in_cache:
            self._save_to_cache(cur=file_key)

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
            return _json_loader(local_file.as_posix())
        return None


class InvalidCorpusError(Exception):
    pass


def default_assembly(stmts):
    from indra.belief.wm_scorer import get_eidos_scorer
    from indra.preassembler.hierarchy_manager import get_wm_hierarchies
    hm = get_wm_hierarchies()
    scorer = get_eidos_scorer()
    stmts = ac.run_preassembly(stmts, belief_scorer=scorer,
                               return_toplevel=True,
                               flatten_evidence=True,
                               normalize_equivalences=True,
                               normalize_opposites=True,
                               normalize_ns='WM',
                               flatten_evidence_collect_from='supported_by',
                               poolsize=4,
                               hierarchies=hm)
    stmts = ac.merge_groundings(stmts)
    stmts = ac.merge_deltas(stmts)
    stmts = ac.standardize_names_groundings(stmts)
    return stmts


def _make_wm_ontology():
    return YamlHierarchyManager(load_yaml_from_url(wm_ont_url),
                                rdf_graph_from_yaml, True)


def _clean_key(s3key):
    # Check if default_key_base ('indra_models') is present in key
    s3key = s3key if default_key_base in s3key else \
        default_key_base + '/' + s3key

    # Replace double slashes
    s3key = s3key.replace('//', '/')

    # Ommit file part of key, assume it ends with json if it is present
    s3key = '/'.join([s for s in s3key.split('/')[:-1]]) if \
        s3key.endswith('.json') else s3key

    # Ensure last char in string is not '/'
    s3key = s3key[:-1] if s3key.endswith('/') else s3key

    return s3key


def _stmts_dict_to_json(stmt_dict):
    """Make a json representation from dict of statements

    This function is the inverse of _json_to_stmts_dict()

    Parameters
    ----------
    stmt_dict : dict
        Dict with statements keyed by their uuid's: {uuid: stmt}

    Returns
    -------
    list(json)
        A list of json statements
    """
    return [s.to_json() for _, s in stmt_dict.items()]


def _json_to_stmts_dict(stmt_jsons):
    """Return dict of statements keyed by uuid's from json statements

    This function is the inverse of _stmts_dict_to_json()

    Parameters
    ----------
    stmt_jsons : list(json)
        A list of json statements

    Returns
    -------
    dict
        Dict with statements keyed by their uuid's: {uuid: stmt}
    """
    loaded_stmts = [Statement._from_json(s) for s in stmt_jsons]
    return {s.uuid: s for s in loaded_stmts}


class LiveCurator(object):
    """Class coordinating the real-time curation of a corpus of Statements.

    Parameters
    ----------
    scorer : indra.belief.BeliefScorer
        A scorer object to use for the curation
    corpora : dict[str, Corpus]
        A dictionary mapping corpus IDs to Corpus objects.
    """

    def __init__(self, scorer=None, corpora=None):
        self.corpora = corpora if corpora else {}
        self.scorer = scorer if scorer else get_eidos_bayesian_scorer()
        self.ont_manager = _make_wm_ontology()
        self.eidos_reader = EidosReader()

    # TODO: generalize this to other kinds of scorers
    def reset_scorer(self):
        """Reset the scorer used for couration."""
        self.scorer = get_eidos_bayesian_scorer()
        for corpus_id, corpus in self.corpora.items():
            corpus.curations = {}

    def get_corpus(self, corpus_id, check_s3=False, use_cache=True):
        """Return a corpus given an ID.

        If the corpus ID cannot be found, an InvalidCorpusError is raised.

        Parameters
        ----------
        corpus_id : str
            The ID of the corpus to return.
        check_s3 : bool
            If True, look on S3 for the corpus if it's not currently loaded.
            Default: False.
        use_cache : bool
            If True, look in local cache before trying to find corpus on s3.
            If True while check_s3 if False, this option will be ignored.
            Default: False.

        Returns
        -------
        Corpus
            The corpus with the given ID.
        """
        try:
            corpus = self.corpora.get(corpus_id)
            if check_s3 and corpus is None:
                logger.info('Corpus not loaded, looking on S3')
                corpus = Corpus.load_from_s3(s3key=corpus_id,
                                             force_s3_reload=not use_cache,
                                             raise_exc=True)
                logger.info('Adding corpus to loaded corpora')
                self.corpora[corpus_id] = corpus
            elif corpus is None:
                raise InvalidCorpusError
            return corpus
        except KeyError:
            raise InvalidCorpusError

    def submit_curation(self, corpus_id, curations):
        """Submit correct/incorrect curations fo a given corpus.

        Parameters
        ----------
        corpus_id : str
            The ID of the corpus to which the curations apply.
        curations : dict
            A dict of curations with keys corresponding to Statement UUIDs and
            values corresponding to correct/incorrect feedback.
        """
        corpus = self.get_corpus(corpus_id, check_s3=True, use_cache=True)
        # Start tabulating the curation counts
        prior_counts = {}
        subtype_counts = {}
        # Take each curation from the input
        for uuid, correct in curations.items():
            # Save the curation in the corpus
            # TODO: handle already existing curation
            stmt = corpus.statements.get(uuid)
            if stmt is None:
                logger.warning('%s is not in the corpus.' % uuid)
                continue
            corpus.curations[uuid] = correct
            # Now take all the evidences of the statement and assume that
            # they follow the correctness of the curation and contribute to
            # counts for their sources
            for ev in stmt.evidence:
                # Make the index in the curation count list
                idx = 0 if correct else 1
                extraction_rule = ev.annotations.get('found_by')
                # If there is no extraction rule then we just score the source
                if not extraction_rule:
                    try:
                        prior_counts[ev.source_api][idx] += 1
                    except KeyError:
                        prior_counts[ev.source_api] = [0, 0]
                        prior_counts[ev.source_api][idx] += 1
                # Otherwise we score the specific extraction rule
                else:
                    try:
                        subtype_counts[ev.source_api][extraction_rule][idx] \
                            += 1
                    except KeyError:
                        if ev.source_api not in subtype_counts:
                            subtype_counts[ev.source_api] = {}
                        subtype_counts[ev.source_api][extraction_rule] = [0, 0]
                        subtype_counts[ev.source_api][extraction_rule][idx] \
                            += 1
        # Finally, we update the scorer with the new curation counts
        self.scorer.update_counts(prior_counts, subtype_counts)

    def save_curations(self, corpus_id, save_to_cache=True):
        """Save the current state of curations for a corpus given its ID

        If the corpus ID cannot be found, an InvalidCorpusError is raised.

        Parameters
        ----------
        corpus_id : str
            the ID of the corpus to save the
        save_to_cache : bool
            If True, also save the current curation to the local cache.
            Default: True.
        """
        # Do NOT use cache or S3 when getting the corpus, otherwise it will
        # overwrite the current corpus
        corpus = self.get_corpus(corpus_id, check_s3=False, use_cache=False)
        corpus.upload_curations(corpus_id, save_to_cache=save_to_cache)

    def update_beliefs(self, corpus_id):
        """Return updated belief scores for a given corpus.

        Parameters
        ----------
        corpus_id : str
            The ID of the corpus for which beliefs are to be updated.

        Returns
        -------
        dict
            A dictionary of belief scores with keys corresponding to Statement
            UUIDs and values to new belief scores.
        """
        # TODO check which options are appropriate for get_corpus
        corpus = self.get_corpus(corpus_id)
        be = BeliefEngine(self.scorer)
        stmts = list(corpus.statements.values())
        be.set_prior_probs(stmts)
        # Here we set beliefs based on actual curation
        for uuid, correct in corpus.curations.items():
            stmt = corpus.statements.get(uuid)
            if stmt is None:
                logger.warning('%s is not in the corpus.' % uuid)
                continue
            stmt.belief = correct
        belief_dict = {st.uuid: st.belief for st in stmts}
        return belief_dict

    def update_groundings(self, corpus_id):
        # TODO check which options are appropriate for get_corpus
        corpus = self.get_corpus(corpus_id)

        # Send the latest ontology and list of concept texts to Eidos
        yaml_str = yaml.dump(self.ont_manager.yaml_root)
        concepts = []
        for stmt in corpus.raw_statements:
            for concept in stmt.agent_list():
                concept_txt = concept.db_refs.get('TEXT')
                concepts.append(concept_txt)
        groundings = self.eidos_reader.reground_texts(concepts, yaml_str)
        # Update the corpus with new groundings
        idx = 0
        for stmt in corpus.raw_statements:
            for concept in stmt.agent_list():
                concept.db_refs['UN'] = groundings[idx]
                idx += 1
        assembled_statements = default_assembly(corpus.raw_statements)
        corpus.statements = {s.uuid: s for s in assembled_statements}
        return assembled_statements


# From here on, a Flask app built around a LiveCurator is implemented

curator = LiveCurator(corpora=corpora)


@app.route('/reset_curation', methods=['POST'])
def reset_curation():
    """Reset the curations submitted until now."""
    if request.json is None:
        abort(Response('Missing application/json header.', 415))
    curator.reset_scorer()
    return jsonify({})


@app.route('/submit_curation', methods=['POST'])
def submit_curation():
    """Submit curations for a given corpus.

    The submitted curations are handled to update the probability model but
    there is no return value here. The update_belief function can be called
    separately to calculate update belief scores.

    Parameters
    ----------
    corpus_id : str
        The ID of the corpus for which the curation is submitted.
    curations : dict
        A set of curations where each key is a Statement UUID in the given
        corpus and each key is 0 or 1 with 0 corresponding to incorrect and
        1 corresponding to correct.
    """
    if request.json is None:
        abort(Response('Missing application/json header.', 415))
    # Get input parameters
    corpus_id = request.json.get('corpus_id')
    curations = request.json.get('curations', {})
    try:
        curator.submit_curation(corpus_id, curations)
    except InvalidCorpusError:
        abort(Response('The corpus_id "%s" is unknown.' % corpus_id, 400))
        return
    return jsonify({})


@app.route('/update_beliefs', methods=['POST'])
def update_beliefs():
    """Return updated beliefs based on current probability model."""
    if request.json is None:
        abort(Response('Missing application/json header.', 415))
    # Get input parameters
    corpus_id = request.json.get('corpus_id')
    try:
        belief_dict = curator.update_beliefs(corpus_id)
    except InvalidCorpusError:
        abort(Response('The corpus_id "%s" is unknown.' % corpus_id, 400))
        return
    return jsonify(belief_dict)


@app.route('/add_ontology_entry', methods=['POST'])
def add_ontology_entry():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    # Get input parameters
    entry = request.json.get('entry')
    examples = request.json.get('examples', [])
    # Add the entry and examples to the in-memory representation
    # of the onotology
    curator.ont_manager.add_entry(entry, examples)
    return jsonify({})


@app.route('/reset_ontology', methods=['POST'])
def reset_ontology():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    # Reload the original ontology
    curator.ont_manager = _make_wm_ontology()

    return jsonify({})


@app.route('/update_groundings', methods=['POST'])
def update_groundings():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    # Get input parameters
    corpus_id = request.json.get('corpus_id')
    # Run the actual regrounding
    stmts = curator.update_groundings(corpus_id)
    stmts_json = stmts_to_json(stmts)
    return jsonify(stmts_json)


if __name__ == '__main__':
    # Process arguments
    parser = argparse.ArgumentParser(
        description='Choose a corpus for live curation.')
    parser.add_argument('--json')
    parser.add_argument('--raw_json')
    parser.add_argument('--pickle')
    parser.add_argument('--corpus_id', default='1')
    parser.add_argument('--host', default='0.0.0.0')
    parser.add_argument('--port', default=8001, type=int)
    parser.add_argument('--aws-cred', type=str, default='default',
                        help='The name of the credential set to use when '
                             'connecting to AWS services. If the name is not '
                             'found in your AWS config, `[default]`  is used.')
    args = parser.parse_args()

    # Load the corpus; If no corpus is provided, raise ValueError
    if args.corpus_id == '1' and (not args.pickle and not args.json):
        raise ValueError('Must specify --corpus_id OR (--pickle or --json)')
    if args.corpus_id:
        curator.corpora[args.corpus_id] = Corpus.load_from_s3(
            s3key=args.corpus_id,
            aws_name=args.aws_cred
        )
        logger.info('Loaded corpus %s from S3 with %d statements and %d '
                    'curation entries' %
                    (args.corpus_id,
                     len(curator.corpora[args.corpus_id].statements),
                     len(curator.corpora[args.corpus_id].curations)))
    else:
        if args.json:
            stmts = stmts_from_json_file(args.json)
        elif args.pickle:
            with open(args.pickle, 'rb') as fh:
                stmts = pickle.load(fh)
        if args.raw_json:
            raw_stmts = stmts_from_json_file(args.raw_json)
        else:
            raw_stmts = None
        logger.info('Loaded corpus from provided file with %d statements.' %
                    len(stmts))
        # If loaded from file, the key will be '1'
        curator.corpora[args.corpus_id] = Corpus(stmts, raw_stmts,
                                                 args.aws_cred)

    # Run the app
    app.run(host=args.host, port=args.port, threaded=False)
