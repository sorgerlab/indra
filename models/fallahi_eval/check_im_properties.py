import pickle
import itertools
from pysb import kappa
from indra.explanation.model_checker import remove_im_params
from assemble_models import prefixed_pkl

print("Loading")
with open(prefixed_pkl('pysb_model'), 'rb') as f:
    model = pickle.load(f)

print("Generating influence map")
im = kappa.influence_map(model)
remove_im_params(model, im)

predecessors = im.predecessors_iter
successors = im.successors_iter
combos = list(itertools.combinations(im.nodes(), 2))
# Examine all pairs of nodes
for ix, (p1, p2) in enumerate(combos):
    if ix % 1000 == 0:
        print("%d of %d" % (ix+1, len(combos)))
    p1_children = set(successors(p1))
    p2_children = set(successors(p2))
    if p1_children == p2_children:
        #print("Same children")
        pass
    elif not p1_children.intersection(p2_children):
        #print("No shared children")
        pass
    # Some shared children, not all
    else:
        #print("Non-equal children sets")
        # Here, the only difference between the sets is that they have mutual
        # edges to each other
        if p1_children.difference(p2_children) == set([p2]) and \
           p2_children.difference(p1_children) == set([p1]):
            pass
            #print("    Mutual activation")
        # This means that the two rules may share a common observable--which is
        # fine, as long as they don't have edges to each other (which would be
        # pruned if they existed)
        # They have a link to each other, but the child
        elif p1 in p2_children or p2 in p1_children:
            print('%s, %s' % (p1, p2))
        # They have partially overlapping sets but not links to each other
        else:
            pass
