import glob
import pickle
import matplotlib.pyplot as plt

fnames = glob.glob('combined/other_outputs/*.pkl')
all_stmts = []
for fn in fnames:
    stmts = pickle.load(open(fn, 'rb'))
    all_stmts += stmts

count = {}
belief = {}
keys = ['reach', 'trips', 'nactem', 'reach_trips', 'nactem_trips',
        'nactem_reach', 'nactem_reach_trips']
for st in all_stmts:
    sources = []
    for ev in st.evidence:
        sources.append(ev.source_api)
    sources = set(sources)
    key = str('_'.join(sorted(list(sources))))
    try:
        count[key] += 1
        if st.belief >= 0.5:
            belief[key].append(st.belief)
    except KeyError:
        count[key] = 1
        if st.belief > 0.5:
            belief[key] = [st.belief]

data = []
for k in keys:
    data.append(belief[k])

plt.boxplot(data)
labels = []
for k in keys:
    k = k.replace('nactem', 'nactem/isi')
    names = k.split('_')
    names = [n.upper() for n in names]
    labels.append('+'.join(names))
labels = ['RE', 'TR', 'NI', 'RE + TR', 'NI + TR', 'NI + RE', 'NI + RE + TR']
plt.ylabel('Belief probability of INDRA Statement')
plt.xticks(range(1,8), labels, rotation='vertical')
plt.show()
