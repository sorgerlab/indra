import csv
from indra.bel import bel_api
from indra.biopax import biopax_api as ba
from indra.preassembler import Preassembler, render_stmt_graph
from indra.preassembler.hierarchy_manager import entity_hierarchy as eh, \
                                                 modification_hierarchy as mh
import pickle
from indra.preassembler.sitemapper import default_mapper as sm
from indra.statements import *
from pysb import kappa
from indra.assemblers import PysbAssembler


class GeneNetwork(object):
    """Build a set of INDRA statements for a given gene list from databases.

    Parameters
    ----------
    gene_list : string
        List of gene names.
    basename : string
        Filename prefix to be used for caching of intermediates (Biopax OWL
        file, pickled statement lists, etc.)
    """

    def __init__(self, gene_list, basename):
        if not gene_list:
            raise ValueError("Gene list must contain at least one element.")
        self.gene_list = gene_list
        self.basename = basename

    def get_bel_stmts(self, filter=False):
        """Get relevant statements from the BEL large corpus.

        Performs a series of neighborhood queries and then takes the union of
        all the statements. Because the query process can take a long time for
        large gene lists, the resulting list of statements are cached in a
        pickle file with the filename `<basename>_bel_stmts.pkl`.  If the
        pickle file is present, it is used by default; if not present, the
        queries are performed and the results are cached.

        Parameters
        ----------
        filter : bool
            If True, includes only those statements that exclusively mention
            genes in :py:attr:`gene_list`. Default is False.

        Returns
        -------
        list of :py:class:`indra.statements.Statement`
            List of INDRA statements extracted from the BEL large corpus.
        """
        bel_stmt_path = '%s_bel_stmts.pkl' % self.basename
        # Check for cached BEL stmt file
        if os.path.isfile(bel_stmt_path):
            print "Loading BEL statements from %s" % bel_stmt_path
            with open(bel_stmt_path) as f:
                bel_statements = pickle.load(f)
        # The file was not found, so run the queries
        else:
            print("No BEL statement file found at %s, running queries" %
                  bel_stmt_path)
            bel_statements = []
            for gene in self.gene_list:
                print "Getting BEL statements for gene", gene
                bel_proc = bel_api.process_ndex_neighborhood([gene])
                bel_statements += bel_proc.statements
            # Save to pickle file
            with open(bel_stmt_path, 'w') as f:
                pickle.dump(bel_statements, f)
        # Optionally filter out statements not involving only our gene set
        if filter:
            print("Filtering statements to match gene list")
            bel_statements = [s for s in bel_statements
                              if all([(agent.name in self.gene_list)
                                      for agent in s.agent_list()])]
        return bel_statements

    def get_biopax_stmts(self, filter=False):
        """Get relevant statements from Pathway Commons.

        Performs a "paths between" query for the genes in :py:attr:`gene_list`
        and uses the results to build statements. This function caches two
        files: the list of statements built from the query, which is cached in
        `<basename>_biopax_stmts.pkl`, and the OWL file returned by the Pathway
        Commons Web API, which is cached in `<basename>_pc_pathsbetween.owl`.
        If these cached files are found, then the results are returned based
        on the cached file and Pathway Commons is not queried again.

        Parameters
        ----------
        filter : bool
            If True, includes only those statements that exclusively mention
            genes in :py:attr:`gene_list`. Default is False.

        Returns
        -------
        list of :py:class:`indra.statements.Statement`
            List of INDRA statements extracted from Pathway Commons.
        """
        # Check for cached Biopax stmt file
        biopax_stmt_path = '%s_biopax_stmts.pkl' % self.basename
        if os.path.isfile(biopax_stmt_path):
            print "Loading Biopax statements from %s" % biopax_stmt_path
            with open(biopax_stmt_path) as f:
                bp_statements = pickle.load(f)
            return bp_statements
        # The statement file was not found, so run the query
        print("No Biopax statement file found at %s, running query" %
              biopax_stmt_path)
        # Check for cached file before querying Pathway Commons Web API
        biopax_ras_owl_path = '%s_pc_pathsbetween.owl' % self.basename
        if os.path.isfile(biopax_ras_owl_path):
            print "Loading Biopax from OWL file", biopax_ras_owl_path
            bp = ba.process_owl(biopax_ras_owl_path)
        # OWL file not found; do query and save to file
        else:
            print "Biopax OWL file not found, querying Pathway Commons web API"
            bp = ba.process_pc_pathsbetween(gene_list)
            bp.save_model(biopax_ras_owl_path)
        # Extract statements from Biopax model
        bp.get_phosphorylation()
        bp.get_dephosphorylation()
        bp.get_acetylation()
        bp.get_palmitoylation()
        bp.get_glycosylation()
        bp.get_activity_modification()
        # Save to pickle file
        with open(biopax_stmt_path, 'w') as f:
            pickle.dump(bp.statements, f)
        # Optionally filter out statements not involving only our gene set
        if filter:
            print("Filtering statements to match gene list")
            bp_statements = [s for s in bp.statements
                              if all([(agent.name in self.gene_list)
                                      for agent in s.agent_list()])]
            return bp_statements
        else:
            return bp.statements

    def get_statements(self, filter=False):
        """Return the combined list of statements from BEL and Pathway Commons.

        Internally calls :py:meth:`get_biopax_stmts` and
        :py:meth:`get_bel_stmts`.

        Parameters
        ----------
        filter : bool
            If True, includes only those statements that exclusively mention
            genes in :py:attr:`gene_list`. Default is False.

        Returns
        -------
        list of :py:class:`indra.statements.Statement`
            List of INDRA statements extracted the BEL large corpus and Pathway
            Commons.
        """
        bp_stmts = self.get_biopax_stmts(filter=filter)
        bel_stmts = self.get_bel_stmts(filter=filter)

        return bp_stmts + bel_stmts

    def run_preassembly(self, stmts, print_summary=True):
        """Run complete preassembly procedure on the given statements.

        Results are returned as a dict and stored in the attribute
        :py:attr:`results`. They are also saved in the pickle file
        `<basename>_results.pkl`.

        Parameters
        ----------
        stmts : list of :py:class:`indra.statements.Statement`
            Statements to preassemble.
        print_summary : bool
            If True (default), prints a summary of the preassembly process to
            the console.

        Returns
        -------
        dict
            A dict containing the following entries:

            - `raw`: the starting set of statements before preassembly.
            - `duplicates1`: statements after initial de-duplication.
            - `valid`: statements found to have valid modification sites.
            - `mapped`: mapped statements (list of
              :py:class:`indra.preassembler.sitemapper.MappedStatement`).
            - `mapped_stmts`: combined list of valid statements and statements
              after mapping.
            - `duplicates2`: statements resulting from de-duplication of the
              statements in `mapped_stmts`.
            - `related2`: top-level statements after combining the statements
              in `duplicates2`.
        """
        # First round of preassembly: remove duplicates before sitemapping
        pa1 = Preassembler(eh, mh, stmts)
        print "Combining duplicates"
        pa1.combine_duplicates()
        # Map sites
        print "Mapping sites"
        (valid, mapped) = sm.map_sites(pa1.unique_stmts)
        # Combine valid and mapped statements into single list
        mapped_stmts = valid + [m.mapped_stmt for m in mapped]
        # Second round of preassembly: de-duplicate and combine related
        pa2 = Preassembler(eh, mh, mapped_stmts)
        print "Combining duplicates again"
        pa2.combine_duplicates()
        pa2.combine_related()
        # Fill out the results dict
        self.results = {}
        self.results['raw'] = stmts
        self.results['duplicates1'] = pa1.unique_stmts
        self.results['valid'] = valid
        self.results['mapped'] = mapped
        self.results['mapped_stmts'] = mapped_stmts
        self.results['duplicates2'] = pa2.unique_stmts
        self.results['related2'] = pa2.related_stmts
        # Print summary
        if print_summary:
            print
            print("Starting number of statements: %d" % len(stmts))
            print("After duplicate removal: %d" % len(pa1.unique_stmts))
            print("Unique statements with valid sites: %d" % len(valid))
            print("Unique statements with invalid sites: %d" % len(mapped))
            print("After post-mapping duplicate removal: %d" %
                  len(pa2.unique_stmts))
            print("After combining related statements: %d" %
                  len(pa2.related_stmts))
        # Save the results
        results_filename = '%s_results.pkl' % self.basename
        with open(results_filename, 'w') as f:
            pickle.dump(self.results, f)
        return self.results


def grounding_filter(stmts):
    """Filter a set of statements to include only those with grounded entities.

    If a statement contains only agents having a non-empty `db_refs` dict,
    it is included in the result

    Parameters
    ----------
    stmts : list of :py:class:`indra.statements.Statement`
        Statements to filter.

    Returns
    -------
    list of :py:class:`indra.statements.Statement`
        Statements after filtering.
    """
    grounded_stmts = []
    for stmt in stmts:
        agents = [a for a in stmt.agent_list() if a is not None]
        if all(a.db_refs for a in agents):
            grounded_stmts.append(stmt)
    return grounded_stmts


if __name__ == '__main__':

    # STEP 0: Get gene list
    gene_list = []
    # Get gene list from ras_pathway_proteins.csv
    with open('../../data/ras_pathway_proteins.csv') as f:
        csvreader = csv.reader(f, delimiter='\t')
        for row in csvreader:
            gene_list.append(row[0].strip())

    gn = GeneNetwork(gene_list, 'ras_genes')
    stmts = gn.get_statements(filter=True)
    grounded_stmts = grounding_filter(stmts)
    results = gn.run_preassembly(grounded_stmts)

    #import cProfile
    #import pstats
    #profile = cProfile.Profile()
    #profile.dump_stats('related_stats')
    #stats = pstats.Stats('related_stats')

"""
sublist = ['RAF1', 'MAP2K1', 'MAPK1', 'KSR1', 'DUSP1', 'KRAS', 'AKT1', 'PDPK1']

subgraph = []
for s in results['related2']:
    in_sub_list = [agent is None or agent.name in sublist
                   for agent in s.agent_list()]
    if all(in_sub_list):
        subgraph.append(s)

with open('ras_genes_subgraph.pkl', 'w') as f:
    pickle.dump(subgraph, f)

# STEP 5 Kappa contact map
with open('ras_genes_subgraph.pkl') as f:
    stmts = pickle.load(f)
stmts = [s for s in stmts \
         if not (s.enz.name == 'MAPK1' and s.sub.name == 'RAF1')]
stmts = [s for s in stmts \
         if not (s.enz.name == 'PDPK1' and s.sub.name == 'AKT1')]


pya = PysbAssembler(policies='interactions_only')
pya.add_statements(stmts)
pya.make_model()
pya.model.name = 'ras_subgraph'
g = kappa.contact_map(pya.model)
g.draw('ras_genes_subgraph_cm.pdf', prog='dot')

# STEP 6
# Load the results from reading 100 papers
with open('pmc_papers/preassembler.pkl') as f:
    braf_stmts = pickle.load(f)

pa3 = Preassembler(eh, mh, braf_stmts)
pa3.combine_duplicates()

# Filter statements based on whether they were grounded to a gene/protein
braf_grounded = []
for s in pa3.unique_stmts:
    hgnc_grounded = [agent is None or
                     agent.db_refs.get('HGNC', None) is not None or
                     agent.db_refs.get('UP', None) is not None
                     for agent in s.agent_list()]
    if all(hgnc_grounded):
        braf_grounded.append(s)

# Filter statements based on whether the entities are in the McCormick list
braf_filtered = []
for s in braf_grounded:
    in_ras_list = [agent is None or agent.name in gene_list
                   for agent in s.agent_list()]
    if all(in_ras_list):
        braf_filtered.append(s)
"""


