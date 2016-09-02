import pickle
import sys
import csv
import numpy as np
from indra.assemblers.english_assembler import EnglishAssembler

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print "Usage: %s stmt_file" % sys.argv[0]
        sys.exit()

    input_file = sys.argv[1]
    with open(input_file) as f:
        stmts_by_paper = pickle.load(f)

    stmts_by_rule = {}
    for paper, stmts in stmts_by_paper.items():
        for stmt in stmts:
            found_by_rule = stmt.evidence[0].annotations['found_by']
            stmt_list = stmts_by_rule.get(found_by_rule)
            if stmt_list is None:
                stmts_by_rule[found_by_rule] = [stmt]
            else:
                stmt_list.append(stmt)

    with open('reach_stmts_by_rule.pkl', 'w') as f:
        pickle.dump(stmts_by_rule, f)

    frequencies = [(k, len(v)) for k, v in stmts_by_rule.items()]
    frequencies.sort(key=lambda x: x[1], reverse=True)
    with open('reach_rule_frequencies.tsv', 'w') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        csvwriter.writerows(frequencies)

    sample_rows = []
    max_sample_size = 20
    for rule, freq in frequencies:
        stmts = stmts_by_rule[rule]
        if max_sample_size < len(stmts):
            sample_stmts = np.random.choice(stmts,
                                            max_sample_size, replace=False)
        else:
            sample_stmts = stmts
        for stmt in sample_stmts:
            for ag in stmt.agent_list():
                if ag is not None:
                    ag.name = ag.db_refs.get('TEXT')
            is_hypothesis = stmt.evidence[0].epistemics.get('hypothesis', '')
            is_direct = stmt.evidence[0].epistemics.get('direct', '')
            # Get the English assembly of the statement
            eng = EnglishAssembler([stmt])
            eng_sentence = eng.make_model()
            if eng_sentence == '':
                eng_sentence = str(stmt)
            sample_rows.append([eng_sentence, is_hypothesis, '', '', '',
                                stmt.evidence[0].pmid,
                                stmt.evidence[0].text, rule, freq, stmt,
                                is_direct])

    with open('stmts_by_rule_to_curate.tsv', 'w') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        csvwriter.writerows(sample_rows)

        # Sample 100 stmts from each rule and put in a .tsv file and upload to
        # Google docs

