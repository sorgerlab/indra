from pysb import *
from pysb import kappa

Model()

Monomer('Pervanadate', ['b'])
Monomer('DUSP', ['b'])
Monomer('MAPK1', ['b', 'T185'], {'T185': ['u', 'p']})

Rule('Pvd_binds_DUSP',
     Pervanadate(b=None) + DUSP(b=None) >>
     Pervanadate(b=1) % DUSP(b=1),
     Parameter('k1', 1))
Rule('Pvd_binds_DUSP_rev',
     Pervanadate(b=1) % DUSP(b=1) >>
     Pervanadate(b=None) + DUSP(b=None),
     Parameter('k2', 1))
Rule('DUSP_binds_MAPK1_phosT185',
     DUSP(b=None) + MAPK1(b=None, T185='p') >>
     DUSP(b=1) % MAPK1(b=1, T185='p'),
     Parameter('k3', 1))
Rule('DUSP_binds_MAPK1_phosT185_rev',
     DUSP(b=1) % MAPK1(b=1, T185='p') >>
     DUSP(b=None) + MAPK1(b=None, T185='p'),
     Parameter('k4', 1))
Rule('DUSP_dephos_MAPK1_at_T185',
     DUSP(b=1) % MAPK1(b=1, T185='p') >>
     DUSP(b=None) % MAPK1(b=None, T185='u'),
     Parameter('k5', 1))

Observable('MAPK1_pT185', MAPK1(T185='p'))

def remove_im_params(model, im):
    for param in model.parameters:
        im.remove_node(param.name)

im = kappa.influence_map(model)
remove_im_params(model, im)
im.draw('pvd_im.pdf', prog='dot')
