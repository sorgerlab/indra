import json
import pickle
import itertools
from indra.util import write_unicode_csv
from indra.assemblers import PysbAssembler, EnglishAssembler, CyJSAssembler
from indra.explanation.model_checker import ModelChecker
import indra.tools.assemble_corpus as ac
import process_data
import make_stmts_for_checking as make_stmts
from assemble_pysb import set_context, add_observables

def get_path_stmts(results, model, stmts):
    all_path_stmts = []
    for drug, target, polarity, value, found_path, paths, flag in results:
        path_stmts = {}
        for path in paths:
            path_stmts1 = stmts_for_path(path, model, stmts)
            for ps in path_stmts1:
                path_stmts[ps.uuid] = ps
        all_path_stmts.append(path_stmts)
    return all_path_stmts

def stmts_for_path(path, model, stmts):
    path_stmts = []
    for path_rule, sign in path:
        for rule in model.rules:
            if rule.name == path_rule:
                stmt = _stmt_from_rule(model, path_rule, stmts)
                path_stmts.append(stmt)
    return path_stmts

def get_path_genes(all_path_stmts):
    path_genes = []
    for path_stmts in all_path_stmts:
        for stmt in path_stmts.values():
            for agent in stmt.agent_list():
                if agent is not None:
                    path_genes.append(agent.name)
    path_genes = sorted(list(set(path_genes)))
    return path_genes

def _stmt_from_rule(model, rule_name, stmts):
    """Return the INDRA Statement corresponding to a given rule by name."""
    stmt_uuid = None
    for ann in model.annotations:
        if ann.predicate == 'from_indra_statement':
            if ann.subject == rule_name:
                stmt_uuid = ann.object
                break
    if stmt_uuid:
        for stmt in stmts:
            if stmt.uuid == stmt_uuid:
                return stmt

def make_cyjs_network(results, model, stmts):
    path_stmts = get_path_stmts(results, model, stmts)
    path_genes = get_path_genes(path_stmts)
    # Get UUIDs to use as filter
    path_uuids = [list(path.keys()) for path in path_stmts]
    all_path_uuids = []
    for p in path_uuids:
        all_path_uuids += p
    #filtered_stmts = ac.filter_gene_list(stmts, path_genes, 'one')
    filtered_stmts = ac.filter_uuid_list(stmts, all_path_uuids)
    ac.dump_statements(filtered_stmts, 'output/korkut_cyjs_model.pkl')
    ca = CyJSAssembler(filtered_stmts)
    cm = ca.make_model()
    ca.set_CCLE_context(['SKMEL28_SKIN'])
    ca.save_json('output/korkut_model')


def make_english_output(results, model, stmts):
    citations = {}
    citation_count = 1
    for source, target, polarity, value, found_path, paths, flag in results:
        cond = 'How does treatment with %s %s %s?' % \
            (source, 'increase' if polarity == 'positive' else
                     'decrease', target)
        print(cond)
        print('=' * len(cond))
        if paths:
            path = paths[0]
            sentences = []
            for i, (path_rule, sign) in enumerate(path):
                for rule in model.rules:
                    if rule.name == path_rule:
                        stmt = _stmt_from_rule(model, path_rule, stmts)
                        if i == 0:
                            sentences.append('%s is a target of %s.' %
                                            (stmt.agent_list()[0].name, source))

                        # Make citations
                        pmids = [ev.pmid for ev in stmt.evidence if ev.pmid]
                        cit_nums = []
                        for pmid in pmids:
                            cit_num = citations.get(pmid)
                            if cit_num is None:
                                citations[pmid] = citation_count
                                cit_num = citation_count
                                citation_count += 1
                            cit_nums.append(cit_num)
                        if cit_nums:
                            cit_nums = sorted(list(set(cit_nums)))
                            cit_str = ' [%s]' % (','.join([str(c) for c
                                                          in cit_nums]))
                        else:
                            cit_str = ''
                        ea = EnglishAssembler([stmt])
                        sentence = ea.make_model()
                        sentence = sentence[:-1] + cit_str + '.'
                        sentences.append(sentence)
            sentences[-1] = sentences[-1][:-1] + \
                ', which is measured by %s.' % target
            text = ' '.join(sentences)
            print('INDRA\'s hypothesis: ' + text)
        elif found_path:
            print('INDRA determined that there exists an explanation but'
                  ' it is intractable to reconstruct.')
        else:
            print('INDRA couldn\'t find an explanation for this observation.')
        print('\n')
    references = 'References\n==========\n'
    for k, v in sorted(citations.items(), key=lambda x: x[1]):
        references += '[%d] https://www.ncbi.nlm.nih.gov/pubmed/%s\n' % (v, k)
    print(references)

def export_json(results, model, stmts):
    """Export a set of paths in JSON format for visualization."""
    json_dict = {}
    for drug, ab, relation, value, path_found, paths, flag in results:
        if json_dict.get(drug) is None:
            json_dict[drug] = {}
        if json_dict[drug].get(ab) is None:
            json_dict[drug][ab] = {}
        for idx, path in enumerate(paths):
            path_stmts = []
            for rule_name, sign in path[:-1]:
                stmt = _stmt_from_rule(model, rule_name, stmts)
                path_stmts.append(stmt.uuid)
            json_dict[drug][ab][idx] = path_stmts
    return json_dict

if __name__ == '__main__':
    print("Processing data")

    data = process_data.read_data(process_data.data_file)
    data_genes = process_data.get_all_gene_names(data)
    ab_map = process_data.get_antibody_map(data)

    print('Loading data statements.')
    data_stmts, data_values = make_stmts.run(dec_thresh=0.5, inc_thresh=1.5)
    all_data_stmts = [values.values() for values in data_stmts.values()]
    all_data_stmts = itertools.chain.from_iterable(all_data_stmts)
    all_data_stmts = list(itertools.chain.from_iterable(all_data_stmts))

    print('We will check the following drug-ab combinations:\n============')
    for drug, stmtd in data_stmts.items():
        print(drug)
        for ab in stmtd.keys():
            print('-'+ ab)

    agent_obs = list(itertools.chain.from_iterable(ab_map.values()))
    # Here we need to cross-reference the antbody map with the data values
    agent_data = {}
    for drug_name, values in data_values.items():
        agent_data[drug_name] = {}
        for ab_name, value in values.items():
            agents = ab_map[ab_name]
            for agent in agents:
                agent_data[drug_name][agent] = value

    base_stmts = ac.load_statements('output/korkut_model_pysb_before_pa.pkl')
    for st in base_stmts:
        st.uuid = str(st.uuid)

    """
    # Merge the sources of statements
    # stmts = manual_stmts + base_stmts
    stmts = base_stmts
    #stmts = manual_stmts

    # Assemble model
    pa = PysbAssembler()
    pa.add_statements(stmts)
    model = pa.make_model()

    with open('korkut_pysb.pkl', 'wb') as f:
        print("Pickling PySB model")
        pickle.dump(pa.model, f)
    """

    with open('korkut_pysb.pkl', 'rb') as f:
        print("Unpickling PySB model")
        model = pickle.load(f)

    # Some parameters up front
    MAX_PATHS_ONE = 5
    MAX_PATHS_ALL = 5
    MAX_PATH_LENGTH = 6


    # Preprocess and assemble the pysb model
    #model = assemble_pysb(combined_stmts, data_genes, '')
    rerun = True
    if rerun:
        mc = ModelChecker(model, all_data_stmts, agent_obs)
        mc.prune_influence_map()

        # Iterate over each drug/ab statement subset
        results = []
        for drug_name, ab_dict in data_stmts.items():
            agent_values = agent_data[drug_name]
            for ab, stmt_list in ab_dict.items():
                value = data_values[drug_name][ab]
                # For each subset, check statements; if any of them checks
                # out, we're good and can move on to the next group
                print("-- Checking the effect of %s on %s --" % (drug_name, ab))
                relation = 'positive' if value > 1 else 'negative'
                path_found = 0
                paths = []
                for stmt in stmt_list:
                    print("Checking: %s" % stmt)
                    result = \
                        mc.check_statement(stmt,
                                           max_paths=MAX_PATHS_ONE,
                                           max_path_length=MAX_PATH_LENGTH)
                    print(result)

                    if result.path_found:
                        path_found = 1
                        if result.paths:
                            paths += result.paths
                    else:
                        print("No path found")

                    # For efficiency, break out of loop as soon as 5 paths
                    # in total are found
                    if len(paths) >= MAX_PATHS_ALL:
                        break
                if paths:
                    print('===========================')
                    print('Scoring a total of %d paths' % len(paths))
                    scored_result = mc.score_paths(paths, agent_values,
                                                   loss_of_function=True)
                    for res in scored_result:
                        print(res[1])
                        for link in res[0]:
                            print('--->', link[0], link[1])
                    paths = [s[0] for s in scored_result]
                    print('===========================')
                results.append((drug_name, ab, relation, value, path_found,
                                paths, result.result_code))
        with open('pathfinding_results.pkl', 'wb') as fh:
            pickle.dump(results, fh)
    else:
        with open('pathfinding_results.pkl', 'rb') as fh:
            results = pickle.load(fh)

    write_unicode_csv('model_check_results.csv', results)
    path_stmts = get_path_stmts(results, model, base_stmts)
    path_genes = get_path_genes(path_stmts)
    #make_english_output(results, model, base_stmts)
    make_cyjs_network(results, model, base_stmts)
    paths_json = export_json(results, model, base_stmts)
