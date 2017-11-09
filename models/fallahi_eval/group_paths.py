import pickle
from indra.databases import hgnc_client
from indra.explanation.model_checker import stmt_from_rule
from indra.databases.context_client import get_protein_expression
from util import prefixed_pkl
from process_data import cell_lines


def group_scored_paths(scored_paths, model, stmts):
    def gene_from_rule(rule_name, agent_ix):
        rule_stmt = stmt_from_rule(rule_name, model, stmts)
        if not rule_stmt:
            print("Could not get stmt for rule %s" % rule_name)
            return None
        agent = rule_stmt.agent_list()[agent_ix]
        gene_id = agent.db_refs.get('HGNC')
        gene_name = hgnc_client.get_hgnc_name(gene_id)
        return gene_name

    path_details = {}
    groups = set()
    for path, score in scored_paths:
        grouped_path = []
        for ix, (rule_name, polarity) in enumerate(path):
            # The last rule should be an observable, which won't have a
            # statement associated with it (for now)
            if ix == len(path) - 1:
                rule_name = path[-2][0]
                gene_name = gene_from_rule(rule_name, 1)
                polarity = 1
            else:
                gene_name = gene_from_rule(rule_name, 0)
            if gene_name:
                grouped_path.append((gene_name, polarity))
        grouped_path = tuple(grouped_path)
        groups.add(grouped_path)
        if grouped_path in path_details:
            path_details[grouped_path].add((tuple(path), score))
        else:
            path_details[grouped_path] = set([(tuple(path), score)])

    print("Flattened %d paths into %d path groups" %
            (len(scored_paths), len(groups)))
    print("Groups: %s" % groups)
    return groups, path_details


def num_paths_by_group(path_details):
    num_paths = {}
    for group, paths in path_details.items():
        if group in num_paths:
            num_paths[group] += len(paths)
        else:
            num_paths[group] = len(paths)
    return num_paths


def print_list(mylist):
    print('\n'.join([str(i) for i in mylist]))


def print_dict(mydict):
    print('\n'.join(['%s: %s' % (k, v) for k, v in mydict.items()]))


def rank_paths(groups, protein_data):
    def rank(path, data):
        score = 0
        # Use set to eliminate duplicate genes which elevates scores
        genes = set([gene for gene, pol in path])
        for gene in genes:
            amt = data[gene] if data[gene] else 0
            score += amt
        return score

    ranked_dict = {}
    for cell_line, data in protein_data.items():
        if data is None:
            continue
        ranked_paths = []
        for path_group in groups:
            score = rank(path_group, data)
            ranked_paths.append((path_group, score))
        ranked_paths.sort(key=lambda x: x[1], reverse=True)
        ranked_dict[cell_line] = ranked_paths
    return ranked_dict


def all_scores(path_details):
    all_scores = {}
    for group in path_details.keys():
        scores = [score for path, score in path_details[group]]
        all_scores[group] = scores
    return all_scores


def top_scores(path_details):
    top_scores = {}
    for group in path_details.keys():
        scores = [score for path, score in path_details[group]]
        top_scores[group] = max(scores)
    ts_tuples = sorted([(k, v) for k, v in top_scores.items()],
                       key=lambda x: (x[1], -len(x[0])), reverse=True)
    return ts_tuples


def print_top_group_scores(scored_paths, model, stmts):
    groups, path_details = group_scored_paths(scored_paths, model, stmts)
    ts = top_scores(path_details)
    last_sign = 1
    for path, score in ts:
        path_str = ''
        for ix, (node, sign) in enumerate(path):
            if ix == 0:
                path_str += node
            else:
                if sign == last_sign:
                    path_str += ' -> %s' % node
                else:
                    path_str += ' -| %s' % node
            last_sign = sign
        print('%s : score %s' % (path_str, score))


if __name__ == '__main__':
    # Run run_task1.py before running this one
    with open(prefixed_pkl('pysb_stmts'), 'rb') as f:
        stmts = pickle.load(f)
    with open('scored_paths.pkl', 'rb') as f:
        (scored_paths, model) = pickle.load(f)

    all_groups = set()
    all_path_details = {}
    for cell_line, drug_dict in scored_paths.items():
        for drug, paths in drug_dict.items():
            groups, path_details = group_scored_paths(paths, model, stmts)
            for pg, path_list in path_details.items():
                if pg in all_path_details:
                    all_path_details[pg] |= path_list
                else:
                    all_path_details[pg] = path_list
            all_groups |= groups

    gene_names = set([tup[0] for group in all_groups for tup in group])
    cell_lines_skin = ['%s_SKIN' % cl for cl in cell_lines]
    protein_data = get_protein_expression(gene_names, cell_lines_skin)

    top_scores = {}
    for group in all_groups:
        scores = [score for path, score in all_path_details[group]]
        top_scores[group] = max(scores)

    print_dict(top_scores)
