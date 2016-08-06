from indra import trips, reach
from indra.literature import id_lookup
from assembly_eval import have_file, run_assembly
from rasmodel import ras_model_statements

if __name__ == '__main__':
    pmc_ids = [s.strip() for s in open('pmcids.txt', 'rt').readlines()]
    pmids = [id_lookup(pmcid)['pmid'] for pmcid in pmc_ids]

    for pmid, pmcid in zip(pmids, pmc_ids):
        print 'Processing %s...' % pmcid
        trips_fname = 'trips/' + pmcid + '.ekb'
        tp = trips.process_xml(open(trips_fname).read())
        for s in tp.statements:
            for e in s.evidence:
                e.pmid = pmid
        reach_fname = 'reach/' + pmcid + '.json'
        rp = reach.process_json_file(reach_fname, citation=pmid)
        rasmodel_stmts = ras_model_statements()
        all_statements = tp.statements + rp.statements + rasmodel_stmts
        run_assembly(all_statements, 'combined', pmcid)
