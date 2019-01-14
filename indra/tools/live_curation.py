import sys
from flask import Flask, request, abort, Response
from indra.belief import wm_scorer, BeliefEngine

scorer = wm_scorer.get_eidos_counts()


app = Flask(__name__)


corpora = {
    '1': {}
    }


@app.route('/update_beliefs', methods=['POST'])
def update_beliefs():
    corpus_id = request.json.get('corpus_id')
    curations = request.json.get('curations')
    return_beliefs = request.json.get('return_beliefs', False)
    prior_counts = {}
    subtype_counts = {}
    corpus = corpora.get(corpus_id)
    for uuid, correct in curations:
        stmt = corpus.get(uuid)
        for ev in stmt.evidence:
            extraction_rule = ev.epistemics.get('found_by')
            if not extraction_rule:
                try:
                    prior_counts[ev.source_api][correct] += 1
                except KeyError:
                    prior_counts[ev.source_api] = [0, 0]
                    prior_counts[ev.source_api][correct] += 1
            else:
                try:
                    prior_counts[ev.source_api][extraction_rule][correct] += 1
                except KeyError:
                    prior_counts[ev.source_api][extraction_rule] = [0, 0]
                    prior_counts[ev.source_api][extraction_rule][correct] += 1

    scorer.update_counts(prior_counts, subtype_counts)
    if not return_beliefs:
        return Response({})
    else:
        be = BeliefEngine(scorer)
        stmts = list(corpus.items())
        be.set_prior_probs(stmts)
        return Response(_get_belief_dict(stmts))


def _get_belief_dict(stmts):
    return {st.uuid: st.belief for st in stmts}


if __name__ == '__main__':
    corpus_path = sys.argv[1]
    with open(corpus_path, 'rb') as fh:
        corpus['1'] = pickle.load(fh)
    app.run()
