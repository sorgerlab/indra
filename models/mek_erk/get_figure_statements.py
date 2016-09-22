from indra.tools.gene_network import GeneNetwork

gene_list = ['MAP2K1', 'MAP2K2', 'MAPK1', 'MAPK3']

gn = GeneNetwork(gene_list, 'mek_erk')
stmts = gn.get_statements()


bp_stmts = gn.get_biopax_stmts(filter=False)
res_bp = gn.run_preassembly(bp_stmts)
uniq_bp = res_bp['duplicates1']
bp_stmts_txt = '\n'.join(sorted(['%s x %d' % (str(s), len(s.evidence))
                 for s in uniq_bp]))
with open('bp_stmts.txt', 'w') as f:
    f.write(bp_stmts_txt)

bel_stmts = gn.get_bel_stmts(filter=True)
res_bel = gn.run_preassembly(bel_stmts)
uniq_bel = res_bel['duplicates1']
bel_stmts_txt = '\n'.join(sorted(['%s x %d' % (str(s), len(s.evidence))
                 for s in uniq_bel]))
with open('bel_stmts.txt', 'w') as f:
    f.write(bel_stmts_txt)




