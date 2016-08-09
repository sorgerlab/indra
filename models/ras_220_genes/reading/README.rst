Description of files in this directory
======================================

get_ras_pmids.py
----------------

Python script that searches Pubmed for a set of PMIDs relevant to the Ras 227
gene set. Also contains a few functions for plotting statistics on the PMIDs
returned by the different queries.

Requires the files:

* ../../data/ras_pathway_proteins.csv: McCormick's Ras gene list.

* pmids_oa_txt.txt : list of PMIDs with PMC Open Access content in plain text
  format.

* pmids_oa_xml.txt : list of PMIDs with PMC Open Access content in XML format.

* pmids_auth_xml.txt : list of PMIDs with Author's manuscript text in XML format.

* pmid_pmcid_doi_map.txt : ???

Produces the following files:

* pmids.pkl: OrderedDict, mapping gene names to a set of PMIDs returned by
  querying Pubmed with the HGNC gene name.

* pmids_from_gene.pkl: OrderedDict, mapping gene names to a set of PMIDs
  returned by querying Entrez Gene for a set of references for the gene.

pmid_pmcid_doi_map.txt
----------------------

Tab-separated variable file with columns: PMID, PMCID, DOI, drawn from the
`PMC-ids.csv` file downloaded from PMC.

