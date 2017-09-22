#{"drug_time_conc_line" : {
#        "meta" : ["StartGene", "EndGene", "StartGene > ... >  EndGene", Lenght, Score],
#            "path" : ["uuid_1", uuids, "uuid_n"]
#                }
#                }
#

def get_paths(scored_paths, model, stmts):
    conc = 0.1
    time = 10
    paths = {}
    for cell_line in scored_paths.keys():
        for drug in scored_paths[cell_line].keys():
            scpaths = scored_paths[cell_line][drug]
            path, score = scpaths[0]
            label = '%s_%s_%s_%s' % (drug, time, conc, cell_line)
            paths[label] = {'meta': [], 'path': []}
            path_stmts = stmts_for_path(path, model, stmts)
            uuids = [stmt.uuid for stmt in path_stmts]
            paths[label]['path'] = uuids
    return paths


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
                assert stmt is not None
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
