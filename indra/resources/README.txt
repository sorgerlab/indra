This folder contains non-python resource files used by various INDRA components.

Files that need to be periodically updated
==========================================

Files that can be generated automatically
-----------------------------------------
bel_chebi_map.tsv
- A table containing mappings from ChEBI names in the BEL namespace to
  ChEBI IDs used in INDRA
- It is created by cross-referencing
  http://resource.belframework.org/belframework/latest-release/equivalence/chebi-ids.beleq
  with
  http://resource.belframework.org/belframework/latest-release/equivalence/chebi.beleq

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

bel_indra_map.tsv
- A table mapping BEL family names to Bioentities entries


Files that don't need periodical updates
========================================
amino_acids.tsv
- A table of the 20 amino acids with their full, short and letter names

index_card_schema.json
- The JSON schema for index cards produced by the IndexCardAssembler
  and consumed by the IndexCardProcessor
