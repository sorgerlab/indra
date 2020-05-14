"""This REST service allows real-time curation and belief updates for
a corpus of INDRA Statements."""
import pickle
import logging
import argparse
from os import path
from flask import Flask, request, jsonify, abort, Response
from indra.statements import stmts_from_json_file, stmts_to_json
from indra.ontology.world.ontology import WorldOntology, wm_ont_url

from . import InvalidCorpusError
from .corpus import Corpus
from .curator import LiveCurator
from .util import _json_loader


logger = logging.getLogger('live_curation')
app = Flask(__name__)
corpora = {}


curator = LiveCurator(corpora=corpora, ont_manager=WorldOntology(wm_ont_url))


# From here on, a Flask app built around a LiveCurator is implemented
@app.route('/download_curation', methods=['POST'])
def download_curation():
    """Download the curations for the given corpus id"""
    if request.json is None:
        abort(Response('Missing application/json header.', 415))
    # Get corpus id, reader name
    corpus_id = request.json.get('corpus_id')
    reader_name = request.json.get('reader', 'all')
    try:
        curation_data = curator.get_curations(corpus_id=corpus_id,
                                              reader=reader_name)
    except InvalidCorpusError:
        abort(Response('The corpus_id "%s" is unknown.' % corpus_id, 400))
        return
    return jsonify(curation_data)


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
    curator.ont_manager = WorldOntology(wm_ont_url)

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


@app.route('/update_metadata', methods=['POST'])
def update_metadata():
    # if request.json is None:
    #     abort(Response('Missing application/json header.', 415))

    # try:
    #     # Get input parameters
    #     corpus_id = request.json.get('corpus_id')
    #     meta_data = request.json.get('meta_data')
    #     curator.update_metadata(corpus_id, meta_data, save_to_cache=True)
    # except InvalidCorpusError:
    #     abort(Response('The corpus_id "%s" is unknown.' % corpus_id, 400))
    #     return
    logger.info('Update metadata is currently disabled')
    return jsonify({'result': 'endpoint disabled'})


@app.route('/save_curation', methods=['POST'])
def save_curations():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    try:
        # Get input parameters
        corpus_id = request.json.get('corpus_id')
        curator.save_curation(corpus_id, save_to_cache=True)
    except InvalidCorpusError:
        abort(Response('The corpus_id "%s" is unknown.' % corpus_id, 400))
        return
    return jsonify({})


if __name__ == '__main__':
    # Process arguments
    parser = argparse.ArgumentParser(
        description='Choose a corpus for live curation.')
    parser.add_argument('--json')
    parser.add_argument('--raw_json')
    parser.add_argument('--pickle')
    parser.add_argument('--meta-json', help='Meta data json file')
    parser.add_argument('--corpus_id')
    parser.add_argument('--host', default='0.0.0.0')
    parser.add_argument('--eidos-url', default='http://localhost:9000')
    parser.add_argument('--port', default=8001, type=int)
    parser.add_argument('--aws-cred', type=str, default='default',
                        help='The name of the credential set to use when '
                             'connecting to AWS services. If the name is not '
                             'found in your AWS config, `[default]`  is used.')
    args = parser.parse_args()

    curator.eidos_url = args.eidos_url

    # Load corpus from S3 if corpus ID is provided
    if args.corpus_id and not args.json and not args.pickle:
        curator.corpora[args.corpus_id] = Corpus.load_from_s3(
            corpus_id=args.corpus_id,
            aws_name=args.aws_cred
        )
        logger.info('Loaded corpus %s from S3 with %d statements and %d '
                    'curation entries' %
                    (args.corpus_id,
                     len(curator.corpora[args.corpus_id].statements),
                     len(curator.corpora[args.corpus_id].curations)))

    elif args.json or args.pickle:
        if not args.corpus_id:
            raise ValueError('Must provide a corpus id when loading files '
                             'locally')
        if args.json:
            stmts = stmts_from_json_file(args.json)
        elif args.pickle:
            with open(args.pickle, 'rb') as fh:
                stmts = pickle.load(fh)
        else:
            stmts = None

        if args.raw_json:
            raw_stmts = stmts_from_json_file(args.raw_json)
        else:
            raw_stmts = None

        if args.meta_json and path.isfile(args.meta_json):
            meta_json_obj = _json_loader(args.meta_json)
        else:
            meta_json_obj = None

        if stmts:
            logger.info('Loaded corpus from provided file with %d '
                        'statements.' % len(stmts))
            # If loaded from file, the key will be '1'
            curator.corpora[args.corpus_id] = Corpus(stmts, raw_stmts,
                                                     meta_json_obj,
                                                     args.aws_cred)

    # Run the app
    app.run(host=args.host, port=args.port, threaded=False)
