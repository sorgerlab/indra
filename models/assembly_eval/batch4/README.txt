The structure of this folder is as follows:
- run_trips_eval.py - runs the assembly evaluation with the TRIPS parser
- run_reach_eval.py - runs the assembly evaluation with the REACH parser
- run_combined.py - runs the assembly evaluation by combining the results of both parsers
- assembly_eval.py - contains the shared code for the above 3 runs, is not meant to be 
                    used directly
Output folders
- reach/ - this folder contains the output of assembly based on REACH parsing
- trips/ - this folder contains the output of assembly based on TRIPS parsing
- combined/ - this folder contains the output of assembly based on combined reader outputs

Output files
- index_cards/PMC*-*.json - Ranked list of generated index cards
- other_outputs/PMC*_statements.tsv - A tab separated spreadsheet of all INDRA
                                      Statements extracted
- other_outputs/PMC*_related_statements.tsv - A tab separated spreadsheet of 
                                              all top-level INDRA Statements 
                                              after running assembly
- other_outputs/PMC*_graph.pdf - A graph showing hierarchical relationships 
                                 between Statements (duplicates not shown)
