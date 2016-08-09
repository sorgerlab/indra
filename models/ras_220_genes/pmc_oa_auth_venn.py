import texttable
import csv
import pickle
import matplotlib_venn as mv
from matplotlib import pyplot as plt
import plot_formatting as pf

pf.set_fig_params()

pmid_map = {}
with open('pmid_pmcid_doi_map.txt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for row in csvreader:
        doi = None if row[2] == '' else row[2]
        pmid_map[row[0]] = (row[1], doi)

with open('pmids_oa_txt.txt') as f:
    pmid_oa_txt = set([line.strip('\n') for line in f.readlines()])

with open('pmids_oa_xml.txt') as f:
    pmid_oa_xml = set([line.strip('\n') for line in f.readlines()])

with open('pmids_auth_xml.txt') as f:
    pmid_auth_xml = set([line.strip('\n') for line in f.readlines()])

plt.figure(figsize=(5, 5), dpi=150)
res = mv.venn3([set(pmid_map.keys()), pmid_oa_xml, pmid_auth_xml],
            ('In PMC with PMID', 'OA subset', 'Author MS'))

for label in res.subset_labels:
    if label is not None:
        label.set_text('{:,d}'.format(int(label.get_text())))

plt.savefig('pmc_oa_auth_venn.png', dpi=150)
plt.savefig('pmc_oa_auth_venn.pdf')

#with open('pmids_from_gene.pkl') as f:
#    pmids_by_gene = pickle.load(f)
#with open('pmids.pkl') as f:
#    pmids_by_name = pickle.load(f)

# Venn diagram of PMC articles
# Articles in Pubmed not in PMC (common--vast majority of articles in PM).
# Articles in PMC not in PM (more rare, but happens ~8% of the time)
# Articles in PMC without an associated DOI (~10%)
# Articles in Pubmed without DOI (?%)
"""

# Num articles in PMC without DOI
num_without_doi = 0
no_doi_set = set([])
for pmid in pmid_map.keys():
    if pmid_map[pmid][1] is None:
        no_doi_set.add(pmid)

total_pmids = len(pmid_map.keys())
print "total", total_pmids
print "no doi in PMC", len(no_doi_set)
print "pct", len(no_doi_set) / float(total_pmids)

#uniq_pmids_by_gene = set([ref for gene, refs in pmids_by_gene.items()
#                                for ref in refs])
#uniq_pmids_by_gene = set([ref for gene, refs in pmids_by_name.items()
#                                for ref in refs])
uniq_pmids = set(pmid_map.keys())
in_oa_xml = uniq_pmids.intersection(pmid_oa_xml)
in_auth_xml = uniq_pmids.intersection(pmid_auth_xml)
in_both = in_oa_xml.intersection(in_auth_xml)
in_txt = uniq_pmids.intersection(pmid_oa_txt)
in_txt_not_xml = in_txt.difference(in_auth_xml).difference(in_oa_xml)
print "total", len(uniq_pmids_by_gene)
print "in auth xml", len(in_auth_xml)
print "in oa xml", len(in_oa_xml)
print "in both", len(in_both)
print "in txt", len(in_txt)
print "in txt no xml", len(in_txt_not_xml)
print "in oa xml, no doi", len(in_oa_xml.intersection(no_doi_set))

print "Total with xml full text", len(in_oa_xml.union(in_auth_xml))

# Query the MySQL database on EC2 for the papers in PMC with OA, AUTH, and TXT
# Do Venn diagrams for
"""

