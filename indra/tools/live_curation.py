"""This REST service allows real-time curation and belief updates for
a corpus of INDRA Statements."""
import sys
import pickle
import logging
import argparse
from flask import Flask, request, jsonify, abort, Response
from indra.belief import wm_scorer, BeliefEngine
from indra.statements import stmts_from_json_file


logger = logging.getLogger('live_curation')
app = Flask(__name__)

default_priors = {'hume': [13, 7], 'cwms': [13, 7], 'sofia': [13, 7]}
scorer = wm_scorer.get_eidos_bayesian_scorer(default_priors)
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


@app.route('/update_beliefs', methods=['POST'])
def update_beliefs():
    if request.json is None:
        abort(Response('Missing application/json header.', 415))

    # Get input parameters
    corpus_id = request.json.get('corpus_id')
    curations = request.json.get('curations', {})
    return_beliefs = request.json.get('return_beliefs', False)

    # Get the right corpus
    try:
        corpus = corpora[corpus_id]
    except KeyError:
        abort(Response('The corpus_id "%s" is unknown.' % corpus_id, 400))
        return

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
    scorer.update_counts(prior_counts, subtype_counts)
    # If not belief return is needed, we just stop here
    if not return_beliefs:
        return jsonify({})
    # Otherwise we rerun the belief calculation on the corpus with
    # the updated scorer and return a dict of beliefs
    else:
        be = BeliefEngine(scorer)
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
        return jsonify(belief_dict)


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
    corpora[args.corpus_id] = Corpus(stmts)

    # Run the app
    app.run(host=args.host, port=args.port)
