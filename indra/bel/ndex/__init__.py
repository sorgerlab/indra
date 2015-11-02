# Example: python query_to_bel.py "ARAF BRAF RAF1" RAF_neighborhood.bel

from . import ndexClient as nc
from . import ndexUtil as util
import sys

def query_to_belscript(proteins,outfile=None):
    myNdex = nc.Ndex("http://ndexbio.org")
    myNet = myNdex.getNeighborhood('9ea3c170-01ad-11e5-ac0f-000c29cb28fb', proteins)
    
    myWrapper = util.NetworkWrapper(myNet,removeNamespace=['MGI','RGD','PFM','NCM','PFR','NCR'])
    bel_script = myWrapper.writeBELScript(outfile)
    return bel_script

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: python query_to_bel.py "protein1 protein2 ..." outfile'
        sys.exit()
    
    proteins = sys.argv[1]
    outfile = sys.argv[2]
    query_to_belscript(proteins,outfile)
 
