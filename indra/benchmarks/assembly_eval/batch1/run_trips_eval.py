from indra import trips
from indra.literature import id_lookup
from assembly_eval import have_file, run_assembly

if __name__ == '__main__':
    pmc_ids = ['PMC1234335', 'PMC3178447', 'PMC3690480',
               'PMC4345513', 'PMC534114']
    pmids = [id_lookup(pmcid, 'pmcid')['pmid'] for pmcid in pmc_ids]
    # Use the existing EKB extractions.
    for pmid, pmcid in zip(pmids, pmc_ids):
        folder = 'trips'
        prefix = folder + '/' + pmcid
        print 'Processing %s...' % pmcid
        tp = trips.process_xml(open(prefix + '-20160503T1152.ekb').read())
        # PMIDs from TRIPS need to be set here because it propagates
        # the PMCID by default
        for s in tp.statements:
            for e in s.evidence:
                e.pmid = pmid
        run_assembly(tp.statements, folder, pmcid)
