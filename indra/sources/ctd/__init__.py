"""This module implements an API and processor to extract INDRA
Statements from the Comparative Toxicogenomics Database (CTD),
see http://ctdbase.org/. It currently extracts chemical-gene, gene-disease,
and chemical-disease relationships. In particular, it extracts the curated
(not inferred) and directional/causal relationships from these subsets."""

from .api import process_from_web, process_dataframe, process_tsv
