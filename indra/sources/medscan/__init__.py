"""
MedScan is Elsevier's proprietary text-mining system for reading the
biological literature. This INDRA module enables processing output
files (in CSXML format) from the MedScan system into INDRA Statements.
"""
from .api import process_file, process_directory, process_file_sorted_by_pmid, \
                 process_directory_statements_sorted_by_pmid
