import pickle
import boolean2
import matplotlib.pyplot as plt
from indra.util import plot_formatting as pf
from indra.assemblers import SifAssembler

def get_sim_avgs(bn_str, nsim=100, nsteps=20, off=None, on=None):
    if off is None:
        off_nodes = []
    else:
        off_nodes = off
    if on is None:
        on_nodes = []
    else:
        on_nodes = on
    coll = boolean2.util.Collector()
    bn_str = boolean2.modify_states(bn_str, turnon=on, turnoff=off)
    model = boolean2.Model(text=bn_str, mode='async')
    for i in range(nsim):
        model.initialize()
        model.iterate(steps=nsteps)
        coll.collect(states=model.states, nodes=model.nodes)
    avgs = coll.get_averages(normalize=True)
    return avgs


st = pickle.load(open('statements.pkl', 'rb'))
sa = SifAssembler(st)
sa.make_model()
bn_str = sa.print_boolean_net('ras_pathway_bn.txt')
# Condition 1
off = []
on = ['Growth_factor_proteins']
avgs = get_sim_avgs(bn_str, off=off, on=on)
jun_basic_noinh = avgs['JUN']
# Condition 2
off = ['MAP2K1', 'MAP2K2']
on = ['Growth_factor_proteins']
avgs = get_sim_avgs(bn_str, off=off, on=on)
jun_basic_inh = avgs['JUN']

st_ext = pickle.load(open('extension.pkl', 'rb'))
sa = SifAssembler(st + st_ext)
sa.make_model()
bn_str = sa.print_boolean_net('ras_pathway_ext_bn.txt')
# Condition 1
off = []
on = ['Growth_factor_proteins']
avgs = get_sim_avgs(bn_str, off=off, on=on)
jun_ext_noinh = avgs['JUN']
# Condition 2
off = ['MAP2K1', 'MAP2K2']
on = ['Growth_factor_proteins']
avgs = get_sim_avgs(bn_str, off=off, on=on)
jun_ext_inh = avgs['JUN']
# Condition 3
off = ['MAP2K1', 'MAP2K2', 'MAPK8', 'MAPK9']
on = ['Growth_factor_proteins']
avgs = get_sim_avgs(bn_str, off=off, on=on)
jun_ext_inh2 = avgs['JUN']

pf.set_fig_params()
plt.figure(figsize=(4,4), dpi=300)
plt.ion()
plt.plot(jun_basic_noinh, 'r', linewidth=2, label='Basic model, no inhibitor')
plt.plot(jun_basic_inh, 'r--', linewidth=2, label='Basic model, MEK inhibitor')
plt.plot(jun_ext_noinh, 'b', linewidth=2, label='Extended model, no inhibitor')
plt.plot(jun_ext_inh, 'b--', linewidth=2, label='Extended model, MEK inhibitor')
plt.plot(jun_ext_inh2, 'g--', linewidth=2, label='Extended model, MEK+JNK inhibitor')
plt.ylim(-0.01, 1.01)
plt.ylabel('Average JUN activity')
plt.xlabel('Time steps')
plt.legend()
pf.format_axis(plt.gca())
plt.show()
