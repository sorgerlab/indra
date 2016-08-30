import pickle
import boolean2
import matplotlib.pyplot as plt
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


st = pickle.load(open('models/ras_pathway/statements.pkl', 'rb'))
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

st_ext = pickle.load(open('models/ras_pathway/extension.pkl', 'rb'))
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

plt.figure()
plt.ion()
plt.plot(jun_basic_noinh, 'r', linewidth=2)
plt.plot(jun_basic_inh, 'r--', linewidth=2)
plt.plot(jun_ext_noinh, 'b', linewidth=2)
plt.plot(jun_ext_inh, 'b--', linewidth=2)
plt.ylim(-0.01, 1.01)
plt.ylabel('Average value')
plt.xlabel('Time steps')
plt.show()
