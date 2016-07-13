from indra import trips
from indra.literature import id_lookup
from assembly_eval import have_file, run_assembly

if __name__ == '__main__':
    pmc_ids = [s.strip() for s in open('pmcids.txt', 'rt').readlines()]
    pmids = [id_lookup(pmcid)['pmid'] for pmcid in pmc_ids]
    # Use the existing EKB extractions.
    for pmid, pmcid in zip(pmids, pmc_ids):
        folder = 'trips'
        prefix = folder + '/' + pmcid
        print 'Processing %s...' % pmcid
        tp = trips.process_xml(open(prefix + '.ekb').read())
        # PMIDs from TRIPS need to be set here because it propagates
        # the PMCID by default
        for s in tp.statements:
            for e in s.evidence:
                e.pmid = pmid
        run_assembly(tp.statements, folder, pmcid)
