"""MedScan is Elsevier's proprietary text-mining system for reading the
biological literature. This INDRA module enables processing output
files (in CSXML format) from the MedScan system into INDRA Statements."""
from indra.sources.medscan.api import process_file
