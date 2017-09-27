from util import *
from indra.statements import *
from indra.assemblers import IndexCardAssembler

stmts = pklload('pysb_stmts')
ica = IndexCardAssembler()

counter = 1
base_path = os.path.join(based, 'index_cards')
for st in stmts:
    if isinstance(st, Modification) and st.enz is None:
        continue
    pmids = [ev.pmid for ev in st.evidence if ev.pmid is not None]
    if pmids:
        pmids = ','.join(['PMID%s' % pm for pm in list(set(pmids))])
    else:
        pmids = 'N/A'
    ica = IndexCardAssembler([st], pmc_override=pmids)
    ica.make_model()
    if ica.cards:
        ica.save_model(os.path.join(base_path, 'index_card_%d.json' % counter))
        counter += 1
