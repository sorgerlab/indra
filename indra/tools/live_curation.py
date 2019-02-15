"""This REST service allows real-time curation and belief updates for
a corpus of INDRA Statements."""
import sys
import pickle
import logging
import argparse
from flask import Flask, request, jsonify, abort, Response
from indra.belief import BeliefEngine
from indra.belief.wm_scorer import get_eidos_bayesian_scorer
from indra.statements import stmts_from_json_file


logger = logging.getLogger('live_curation')
app = Flask(__name__)
corpora = {}


class Corpus(object):
    """Represent a corpus of statements with curation."""
    def __init__(self, statements):
        self.statements = {st.uuid: st for st in statements}
        self.curations = {}

    def __str__(self):
        return 'Corpus(%s -> %s)' % (str(self.statements), str(self.curations))

    def __repr__(self):
        return str(self)


class InvalidCorpusError(Exception):
    pass


class LiveCurator(object):
    """Class coordinating the real-time curation of a corpus of Statements."""

    default_priors = {'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]}

    def __init__(self, scorer, corpora=None):
        self.scorer = scorer
        self.corpora = corpora if corpora else {}

    def reset_scorer(self):
        self.scorer = wm_scorer.get_eidos_bayesian_scorer(self.default_priors)

    def get_corpus(self, corpus_id):
        try:
            corpus = self.corpora[corpus_id]
            return corpus
        except KeyError:
            raise InvalidCorpusError

    def submit_curation(self, corpus_id, curations):
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
                        subtype_counts[ev.source_api][extraction_rule][idx] += 1
                    except KeyError:
                        if ev.source_api not in subtype_counts:
                            subtype_counts[ev.source_api] = {}
                        subtype_counts[ev.source_api][extraction_rule] = [0, 0]
                        subtype_counts[ev.source_api][extraction_rule][idx] += 1
        # Finally, we update the scorer with the new curation counts
        self.scorer.update_counts(prior_counts, subtype_counts)

    def update_beliefs(self, corpus_id):
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


curator = LiveCurator(
    scorer=get_eidos_bayesian_scorer(LiveCurator.default_priors),
    corpora=corpora)


@app.route('/reset_curation', methods=['POST'])
def reset_curation():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))
    curator.reset_scorer()
    return jsonify({})


@app.route('/submit_curation', methods=['POST'])
def submit_curation():
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


@app.route('/add_ontlogy_entry', methods=['POST'])
def add_ontology_entry():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    # Get input parameters
    entry = request.json.get('entry')
    examples = request.json.get('examples', [])

    # Add the entry and examples to the in-memory representation
    # of the onotology
    return jsonify({})


@app.route('/reset_ontology', methods=['POST'])
def reset_ontology():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    # Reload the original ontology

    return jsonify({})


@app.route('/update_groundings', methods=['POST'])
def update_groundings():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    # Get input parameters
    corpus_id = request.json.get('corpus_id')

    # Send the latest ontology and list of concept texts to Eidos
    concepts = []
    for uuid, stmt in corpora.get(corpus_id).items():
        for concept in stmt.agent_list():
            concept_txt = concept.db_refs['TEXT']
            concepts.append(concept_txt)
    # Update the corpus with new groundings

    return jsonify({})


@app.route('/run_assembly', methods=['POST'])
def run_assembly():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    # Get input parameters
    corpus_id = request.json.get('corpus_id')

    # Run preassembly

    # Return assembled statement corpus
    return jsonify({})


if __name__ == '__main__':
    # Process arguments
    parser = argparse.ArgumentParser(description='Choose a corpus for live curation.')
    parser.add_argument('--json')
    parser.add_argument('--pickle')
    parser.add_argument('--corpus_id', default='1')
    parser.add_argument('--host', default='0.0.0.0')
    parser.add_argument('--port', default=8001, type=int)
    args = parser.parse_args()

    # Load the corpus
    if args.json:
        stmts = stmts_from_json_file(args.json)
    elif args.pickle:
        with open(args.pickle, 'rb') as fh:
            stmts = pickle.load(fh)
    logger.info('Loaded corpus %s with %d statements.' %
                (args.corpus_id, len(stmts)))
    curator.corpora[args.corpus_id] = Corpus(stmts)

    # Run the app
    app.run(host=args.host, port=args.port)
