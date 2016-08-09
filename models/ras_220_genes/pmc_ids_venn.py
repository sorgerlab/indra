"""Plot a Venn diagram showing the IDs associated with articles in PMC."""

import matplotlib_venn as mv
import csv
from matplotlib import pyplot as plt
import plot_formatting as pf

pf.set_fig_params()

all_pmcids = set([])
has_doi = set([])
has_pmid = set([])

with open('PMC-ids.csv') as f:
    csvreader = csv.reader(f, delimiter=',')
    for row in csvreader:
        pmcid = row[8].strip()
        pmid = row[9].strip()
        doi = row[7].strip()
        all_pmcids.add(pmcid)
        if doi:
            has_doi.add(pmcid)
        if pmid:
            has_pmid.add(pmcid)
print len(all_pmcids)

plt.figure(figsize=(4, 4), dpi=150)
res = mv.venn2([has_doi, has_pmid],
                  ("DOI", "PMID"))
plt.title('IDs for articles in PMC')
num_neither = len(all_pmcids.difference(has_doi).difference(has_pmid))

def commafy(text):
    text_with_commas = ''
    for ix, char in enumerate(reversed(str(text))):
        if ix % 3 == 0 and ix != 0:
            text_with_commas += ','
        text_with_commas += char
    return text_with_commas[::-1]

plt.text(-0.55, -0.8, '(plus %s with no DOI or PMID)' % commafy(num_neither))

# Add commas
for label in res.subset_labels:
    text = str(label.get_text())
    label.set_text(commafy(text))

plt.show()
plt.savefig('pmc_ids_venn.png', dpi=150)
plt.savefig('pmc_ids_venn.pdf')

