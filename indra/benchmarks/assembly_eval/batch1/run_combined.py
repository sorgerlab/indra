from indra import trips, reach
from indra.literature import id_lookup
from assembly_eval import have_file, run_assembly

if __name__ == '__main__':
    pmc_ids = ['PMC1234335', 'PMC3178447', 'PMC3690480',
               'PMC4345513', 'PMC534114']
    pmids = [id_lookup(pmcid, 'pmcid')['pmid'] for pmcid in pmc_ids]

    for pmid, pmcid in zip(pmids, pmc_ids):
        print 'Processing %s...' % pmcid
        trips_fname = 'trips/' + pmcid + '-20160503T1152.ekb'
        tp = trips.process_xml(open(trips_fname).read())
        for s in tp.statements:
            for e in s.evidence:
                e.pmid = pmid
        reach_fname = 'reach/' + pmcid + '.json'
        rp = reach.process_json_file(reach_fname)
        all_statements = tp.statements + rp.statements
        run_assembly(all_statements, 'combined', pmcid)
