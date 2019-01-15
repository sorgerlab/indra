import sys
import pickle
from flask import Flask, request, jsonify
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
    for uuid, correct in curations.items():
        stmt = corpus.get(uuid)
        for ev in stmt.evidence:
            print(ev)
            extraction_rule = ev.annotations.get('found_by')
            print(extraction_rule)
            if not extraction_rule:
                try:
                    prior_counts[ev.source_api][correct] += 1
                except KeyError:
                    prior_counts[ev.source_api] = [0, 0]
                    prior_counts[ev.source_api][correct] += 1
            else:
                try:
                    subtype_counts[ev.source_api][extraction_rule][correct] += 1
                except KeyError:
                    if ev.source_api not in subtype_counts:
                        subtype_counts[ev.source_api] = {}
                    subtype_counts[ev.source_api][extraction_rule] = [0, 0]
                    subtype_counts[ev.source_api][extraction_rule][correct] += 1
    print(scorer.subtype_counts)
    scorer.update_counts(prior_counts, subtype_counts)
    print(scorer.subtype_counts)
    if not return_beliefs:
        return jsonify({})
    else:
        be = BeliefEngine(scorer)
        stmts = list(corpus.values())
        be.set_prior_probs(stmts)
        belief_dict = _get_belief_dict(stmts)
        print(belief_dict)
        return jsonify(belief_dict)


def _get_belief_dict(stmts):
    return {st.uuid: st.belief for st in stmts}


if __name__ == '__main__':
    corpus_path = sys.argv[1]
    with open(corpus_path, 'rb') as fh:
        stmts = pickle.load(fh)
        corpora['1'] = {st.uuid: st for st in stmts}
    print(corpora)
    app.run()
