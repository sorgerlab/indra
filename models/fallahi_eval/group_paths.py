import pickle
from indra.databases import hgnc_client
from util import prefixed_pkl

def group_scored_paths(scored_paths, model, stmts):
    def gene_from_rule(rule_name, agent_ix):
        rule_stmt = _stmt_from_rule(model, rule_name, stmts)
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
                rule_name = path[-2]
                gene_name = gene_from_rule(rule_name, 1)
                polarity = 1
            else:
                gene_name = gene_from_rule(rule_name, 0)
            if gene_name:
                grouped_path.append((gene_name, polarity))
        grouped_path = tuple(grouped_path)
        groups.add(grouped_path)
        if grouped_path in path_details:
            path_details[grouped_path].append((path, score))
        else:
            path_details[grouped_path] = [(path, score)]

    print("Flattened %d paths into %d path groups" %
            (len(scored_paths), len(groups)))
    print("Groups: %s" % groups)
    return groups, path_details

if __name__ == '__main__':
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
                    all_path_details[pg] += path_list
                else:
                    all_path_details[pg] = path_list
            all_groups |= groups


