from indra import trips, reach
from indra.literature import id_lookup
from assembly_eval import have_file, run_assembly
import rasmodel

if __name__ == '__main__':
    pmc_ids = [s.strip() for s in open('pmcids.txt', 'rt').readlines()]
    pmids = [id_lookup(pmcid, 'pmcid')['pmid'] for pmcid in pmc_ids]

    for pmid, pmcid in zip(pmids, pmc_ids):
        print 'Processing %s...' % pmcid
        # Process TRIPS
        trips_fname = 'trips/' + pmcid + '.ekb'
        tp = trips.process_xml(open(trips_fname).read())
        # Force the correct PMID from outside
        for s in tp.statements:
            for e in s.evidence:
                e.pmid = pmid
        # Process REACH
        reach_fname = 'reach/' + pmcid + '.json'
        rp = reach.process_json_file(reach_fname, citation=pmid)
        # Get prior statements
        rasmodel_stmts = rasmodel.get_statements()
        # Combine all statements
        all_statements = tp.statements + rp.statements + rasmodel_stmts
        # Run assembly
        run_assembly(all_statements, 'combined', pmcid)
