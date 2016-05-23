import pickle
from matplotlib import pyplot as plt
import numpy as np
import plot_formatting as pf
from collections import Counter
from indra.statements import *
from indra.preassembler import Preassembler

pf.set_fig_params()

print "Loading REACH results"
with open('reach_stmts.pkl') as f:
    results = pickle.load(f)

counts_per_paper = [(pmid, len(stmts)) for pmid, stmts in results.items()]
zero_pmids = [pmid for pmid, stmts in results.items() if len(stmts) == 0]
counts = np.array([tup[1] for tup in counts_per_paper])

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

plt.ion()

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

