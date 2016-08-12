import pickle
from matplotlib import pyplot as plt
import numpy as np
from indra.tools import plot_formatting as pf
from collections import Counter
from indra.statements import *
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from copy import deepcopy
from indra.databases import uniprot_client
import sys
from indra.preassembler.grounding_mapper import GroundingMapper, \
                                                default_grounding_map
pf.set_fig_params()


if len(sys.argv) < 2:
    print "Usage: %s reach_stmts_file" % sys.argv[0]
    sys.exit()

stmts_file = sys.argv[1]

print "Loading REACH results"
with open(stmts_file) as f:
    results = pickle.load(f)

counts_per_paper = [(pmid, len(stmts)) for pmid, stmts in results.items()]
zero_pmids = [pmid for pmid, stmts in results.items() if len(stmts) == 0]
counts = np.array([tup[1] for tup in counts_per_paper])

# Preassemble to remove duplicates

all_stmts = [stmt for paper_stmts in results.values()
                  for stmt in paper_stmts]

gm = GroundingMapper(default_grounding_map)
map_stmts = gm.map_agents(all_stmts)

ren_stmts = gm.rename_agents(map_stmts)

# Combining duplicates
pa = Preassembler(hierarchies, ren_stmts)
pa.combine_duplicates()

# Sorted by evidence
sorted_stmts = sorted(pa.unique_stmts, key=lambda x: len(x.evidence))


sys.exit()


sys.exit()


"""
with open('grounding_results.csv', 'w') as f:
    for group in grouped_by_text:
        text_string = group[0]
        line = [text_string]
        for db, id, count in group[1]:
            if db == 'UP':
                name = uniprot_client.get_mnemonic(id)
            else:
                name = ''
            line = '%s\t%s\t%s\t%s\t%s\n' % (text_string, db, id, count, name)
            f.write(line)
"""

with open('grounding_human.csv', 'w') as f:
    for group in grouped_by_text:
        text_string = group[0]
        for db, id, count in group[1]:
            if db == 'UP':
                name = uniprot_client.get_mnemonic(id)

            else:
                name = ''
            line = '%s\t%s\t%s\t%s\t%s\n' % (text_string, db, id, count, name)
            f.write(line)



sys.exit()


# Combining duplicates
pa = Preassembler(hierarchies, all_stmts)
pa.combine_duplicates()

# Sorted by evidence
sorted_stmts = sorted(pa.unique_stmts, key=lambda x: len(x.evidence))

#with open('reach_stmts_sorted.pkl') as f:
#    sorted_stmts = pickle.load(f)

# Distribution of pieces of evidence
plt.ion()
plt.figure(figsize=(2, 2), dpi=300)
ax = plt.gca()
plt.plot([len(stmt.evidence) for stmt in sorted_stmts])
pf.format_axis(ax)
ax.set_xlabel('Statement index')
ax.set_ylabel('No. of papers')
ax.set_yscale('log')
plt.subplots_adjust(left=0.23, bottom=0.16)
plt.savefig('reach_evidence_dist.pdf')

# Sorted by unique PMIDs
def uniq_refs_in_evidence(stmt):
    refs = set([])
    for ev in stmt.evidence:
        refs.add(ev.pmid)
    return len(refs)

sorted_uniq_pmids = deepcopy(sorted_stmts)
sorted_uniq_pmids = sorted(sorted_uniq_pmids, key=uniq_refs_in_evidence,
                           reverse=True)

# List of all entities

# Agent counter
agent_counter = Counter(agents)
agent_counter = sorted(agent_counter.items(), key=lambda x: x[1],
                       reverse=True)

# What fraction of statements grounded?
grounded_dict = {}
for paper_stmts in results.values():
    for stmt in paper_stmts:
        stmt_type = type(stmt).__name__
        total_counter = grounded_dict.get((stmt_type, 'total'), 0)
        grounded_dict[(stmt_type, 'total')] = total_counter + 1
        if all([True if (ag is not None and ag.db_refs) else False
                for ag in stmt.agent_list()]):
            grounded_counter = grounded_dict.get((stmt_type, 'grounded'), 0)
            grounded_dict[(stmt_type, 'grounded')] = grounded_counter + 1

phos = [s for paper_stmts in results.values()
          for s in paper_stmts
          if isinstance(s, Phosphorylation)]

# FIG 1: Distribution of numbers of statements
plt.figure(figsize=(2, 2), dpi=300)
plt.hist(counts, bins=20)
ax = plt.gca()
ax.set_xticks([0, 100, 200, 300, 400])
pf.format_axis(ax)
plt.subplots_adjust(left=0.23, bottom=0.16)
ax.set_xlabel('No. of statements')
ax.set_ylabel('No. of papers')
plt.savefig('reach_stmts_per_paper.pdf')

# FIG 2: Number of statements of different types
stmt_types = [type(stmt) for paper_stmts in results.values()
                       for stmt in paper_stmts]
stmt_type_counter = Counter(stmt_types)

fig = plt.figure(figsize=(2, 3), dpi=300)
ax = fig.gca()
sorted_counts = sorted(stmt_type_counter.items(), key=lambda x: x[1],
                       reverse=True)
labels = [t.__name__ for t, n in sorted_counts]
ax.bar(range(len(sorted_counts)), [x[1] for x in sorted_counts])
ax.set_xticks(np.array(ax.get_xticks()) + 0.4)
ax.set_xticklabels(labels, rotation=90)
pf.format_axis(ax)
ax.set_ylabel('No. of statements')
plt.subplots_adjust(left=0.29, bottom=0.34)
plt.savefig('reach_stmt_types.pdf')


