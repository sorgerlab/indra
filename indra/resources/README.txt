This folder contains non-python resource files used by various INDRA components.

Files that need to be periodically updated
==========================================

Files that can be generated automatically
-----------------------------------------
famplex/*.csv
- A set of csv files downloaded from the FamPlex resource at
https://github.com/sorgerlab/famplex

bel_chebi_map.tsv
- A table containing mappings from ChEBI names in the BEL namespace to
  ChEBI IDs used in INDRA
- It is created by cross-referencing
  http://resource.belframework.org/belframework/latest-release/equivalence/chebi-ids.beleq
  with
  http://resource.belframework.org/belframework/latest-release/equivalence/chebi.beleq

famplex_map.tsv
- A table containing mappings from outside name spaces into FamPlex.
- It is generated from famplex/equivalences.csv

cellular_components.tsv
- A table of GO IDs and standardized names corresponding to cellular locations
- It is exported from http://purl.obolibrary.org/obo/go.owl

chebi_to_chembl.tsv
- A table of ChEBI IDs and their corresponding ChEBML IDs
- It is exported from ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/relations.tsv.gz

chebi_to_pubchem.tsv
- A table of ChEBI IDs and their corresponding PUBCHEM IDs
- It is exported from ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/relations.tsv.gz

hgnc_entries.txt
- A table containing HGNC IDs, symbols and mappings to Entrez and UniProt
- Download link: http://tinyurl.com/gnv32vh

kinases.tsv
- A list of human kinases with their Uniprot IDs and HGNC symbols
- Download link: http://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+database%3A%28type%3Ainterpro+IPR000719%29&sort=score

ncit_map.tsv
- A table of NCIT IDs and their corresponding HGNC, CHEBI or GO IDs.
- It is based on https://ncit.nci.nih.gov/ncitbrowser/ajax?action=export_mapping&dictionary=NCIt_to_HGNC_Mapping&version=1.0,
https://ncit.nci.nih.gov/ncitbrowser/ajax?action=export_mapping&dictionary=GO_to_NCIt_Mapping&version=1.1 and
https://ncit.nci.nih.gov/ncitbrowser/ajax?action=export_mapping&dictionary=NCIt_to_ChEBI_Mapping&version=1.0

uniprot_entries.tsv
- A table of all reviewed UniProt entries plus unreviewed human entries
  with UniProt ID, gene name, organism ID and UniProt mnemonic
- Download links: 
    - http://www.uniprot.org/uniprot/?query=reviewed%3Ayes&sort=id&desc=no&columns=id,genes(PREFERRED),organism-id,entry%20name
    - http://www.uniprot.org/uniprot/?sort=id&desc=no&compress=no&query=reviewed%3Ano&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22&force=no&format=tab&columns=id,genes(PREFERRED),organism-id,entry%20name

uniprot_sec_ac.txt
- A table of secondary (deprecated) UniProt accession numbers with mappings to
  primary accession numbers
- Download link: ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/docs/sec_ac.txt

uniprot_subcell_loc.tsv
- A table containing UniProt subcellular locations and their descriptions
- Download link: http://www.uniprot.org/locations/?sort=&desc=&compress=no&query=&force=no&preview=true&format=tab&columns=id

[name]_hierarchy.rdf
- The INDRA hierarchy graphs used for preassembly

Files that are manually curated
-------------------------------
curated_site_map.tsv
- A table containing mappings for amino acid sites on proteins
  that are often misannotated

transcription_factors.tsv
- manually downloaded from http://fantom.gsc.riken.jp/5/sstar/Browse_Transcription_Factors_hg19

grounding_agents.json
- manually curated list of grounding mappings to INDRA Agents with states 
(phosphorylation, mutation, etc.)

cas_to_chebi.tsv
- Manually curated based on common occurrences in Pathway Commons data. Could
be replaced by a more comprehensive map.

phosphatases.tsv
- A list of human phospharases with their HGNC symbols and IDs
- Compiled based on: http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0049943.s009

Files that don't need periodical updates
========================================
amino_acids.tsv
- A table of the 20 amino acids with their full, short and letter names

index_card_schema.json
- The JSON schema for index cards produced by the IndexCardAssembler
  and consumed by the IndexCardProcessor
