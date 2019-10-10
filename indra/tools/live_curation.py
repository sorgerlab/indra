"""This REST service allows real-time curation and belief updates for
a corpus of INDRA Statements."""
import json
import yaml
import boto3
import pickle
import logging
import argparse
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
default_base_name = 'indra_models'
default_profile = 'wm'


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
    def load_from_s3(cls, name, aws_name=default_profile,
                     bucket=default_bucket, key_base_name=default_base_name):
        corpus = cls([], aws_name=aws_name)
        corpus.s3_get(name, bucket, key_base_name)
        return corpus

    def s3_put(self, name, bucket=default_bucket,
               key_base_name=default_base_name):
        """Push a corpus object to S3 in the form of three json files

        The json files representing the object have S3 keys of the format
        <key_base_name>/<name>/<file>.json

        Parameters
        ----------
        name : str
            The name of the model to upload. Is part of the S3 key.
        bucket : str
            The S3 bucket to upload the Corpus to. Default: 'world-modelers'.
        key_base_name : str
            The base object path to upload the json files to. Is part of the
            S3 key. Default: 'indra_models'.

        Returns
        -------
        keys : tuple(str)
            A tuple of three strings giving the S3 key to the pushed objects
        """
        key_base = key_base_name + '/' + name + '/'
        key_base = key_base.replace('//', '/')  # # replace double slashes
        keys = tuple(key_base + s + '.json' for s in ['raw_statements',
                                                      'statements',
                                                      'curations'])
        try:
            s3 = self._get_s3_client()
            # Structure and upload raw statements
            logger.info('Uploading %s to S3' % keys[0])
            s3.put_object(
                Body=json.dumps(stmts_to_json(self.raw_statements)),
                Bucket=bucket, Key=keys[0])

            # Structure and upload assembled statements
            logger.info('Uploading %s to S3' % keys[1])
            s3.put_object(
                Body=_stmts_dict_to_json_str(self.statements),
                Bucket=bucket, Key=keys[1])

            # Structure and upload curations
            logger.info('Uploading %s to S3' % keys[2])
            s3.put_object(
                Body=json.dumps(self.curations),
                Bucket=bucket, Key=keys[2])
            return keys
        except Exception as e:
            logger.exception('Failed to put on s3: %s' % e)
            return None

    def s3_get(self, name, bucket=default_bucket,
               key_base_name=default_profile):
        """Fetch a corpus object from S3 in the form of three json files

        The json files representing the object have S3 keys of the format
        <key_base_name>/<name>/<file>.json

        Parameters
        ----------
        name : str
            The name of the model to fetch. Is part of the S3 key.
        bucket : str
            The S3 bucket to fetch the Corpus from. Default: 'world-modelers'.
        key_base_name : str
            The base object path to fetch the json files from. Is part of the
            S3 key. Default: 'indra_models'.

        """
        key_base = key_base_name + '/' + name + '/'
        key_base = key_base.replace('//', '/')  # replace double slashes
        try:
            s3 = self._get_s3_client()
            # Get and process raw statements
            raw_stmt_jsons = s3.get_object(
                Bucket=bucket,
                Key=key_base+'raw_statements.json')['Body'].read()
            self.raw_statements = stmts_from_json(json.loads(raw_stmt_jsons))

            # Get and process assembled statements from list to dict
            self.statements = _json_str_to_stmts_dict(s3.get_object(
                Bucket=bucket,
                Key=key_base+'statements.json')['Body'].read())

            # Get and process curations if any
            curation_jsons = json.loads(s3.get_object(
                Bucket=bucket,
                Key=key_base+'curations.json')['Body'].read())
            self.curations = {uid: c for uid, c in curation_jsons.items()}

        except Exception as e:
            logger.exception('Failed to get from s3: %s' % e)


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


def _stmts_dict_to_json_str(stmt_dict):
    """Make a json representation from dict of statements keyed by their uuid's

    This function is the inverse of _json_str_to_stmts_dict()

    Parameters
    ----------
    stmt_dict : dict
        Dict with statements keyed by their uuid's: {uuid: stmt}

    Returns
    -------
    json_str : str
        A json compatible string
    """
    return json.dumps([s.to_json() for _, s in stmt_dict.items()])


def _json_str_to_stmts_dict(json_str):
    """Make a dict of statements keyed by their uuid's from json representation

    This function is the inverse of _stmts_dict_to_json_str()

    Parameters
    ----------
    json_str : str
        A json compatible string

    Returns
    -------
    stmt_dict : dict
        Dict with statements keyed by their uuid's: {uuid: stmt}
    """
    stmt_jsons = json.loads(json_str)
    stmts = [Statement._from_json(s) for s in stmt_jsons]
    return {s.uuid: s for s in stmts}


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

    def get_corpus(self, corpus_id):
        """Return a corpus given an ID.

        If the corpus ID cannot be found, an InvalidCorpusError is raised.

        Parameters
        ----------
        corpus_id : str
            The ID of the corpus to return.

        Returns
        -------
        Corpus
            The corpus with the given ID.
        """
        try:
            corpus = self.corpora[corpus_id]
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
        corpus = self.get_corpus(corpus_id)
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

    # Load the corpus
    if args.json:
        stmts = stmts_from_json_file(args.json)
    elif args.pickle:
        with open(args.pickle, 'rb') as fh:
            stmts = pickle.load(fh)
    if args.raw_json:
        raw_stmts = stmts_from_json_file(args.raw_json)
    else:
        raw_stmts = None

    logger.info('Loaded corpus %s with %d statements.' %
                (args.corpus_id, len(stmts)))
    curator.corpora[args.corpus_id] = Corpus(stmts, raw_stmts, args.aws_cred)

    # Run the app
    app.run(host=args.host, port=args.port, threaded=False)
