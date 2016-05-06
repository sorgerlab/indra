from indra import trips
from assembly_eval import have_file, run_assembly

if __name__ == '__main__':
    pmc_ids = ['PMC1234335', 'PMC3178447', 'PMC3690480',
               'PMC4345513', 'PMC534114']

    # Use the existing EKB extractions.
    for pmcid in pmc_ids:
        folder = 'trips'
        prefix = folder + '/' + pmcid
        print 'Processing %s...' % pmcid
        tp = trips.process_xml(open(prefix + '-20160503T1152.ekb').read())

        run_assembly(tp.statements, folder, pmcid)
