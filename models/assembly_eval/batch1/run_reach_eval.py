import shutil
from indra import reach
from indra.literature import pmc_client, get_full_text, id_lookup
from assembly_eval import have_file, run_assembly

if __name__ == '__main__':
    folder = 'reach'
    pmc_ids = ['PMC1234335', 'PMC3178447', 'PMC3690480',
               'PMC4345513', 'PMC534114']
    rerun = False

    # Download the papers if they are not available yet
    pmids = []
    for pmcid in pmc_ids:
        prefix = folder + '/' + pmcid
        if not have_file(prefix + '.nxml') and\
           not have_file(prefix + '.txt'):
            txt, txt_format = get_full_text(pmcid)
            if txt_format == 'nxml':
                fname = prefix + '.nxml'
            else:
                fname = prefix + '.txt'
            with open(fname, 'wt') as fh:
                fh.write(txt.encode('utf-8'))
        pmids.append(id_lookup(pmcid)['pmid'])


    # Read each paper if it hasn't been read yet.
    # Otherwise use the existing json extractions.
    for pmcid, pmid in zip(pmc_ids, pmids):
        prefix = folder + '/' + pmcid
        print 'Processing %s...' % pmcid
        # If REACH already processed it then don't run it again
        if rerun or not have_file(prefix + '.json'):
            if have_file(prefix + '.txt'):
                txt = open(prefix + '.txt').read().decode('utf-8')
                rp = reach.process_text(txt, citation=pmid)
            elif have_file(prefix + '.nxml'):
                rp = reach.process_nxml_file(prefix + '.nxml', citation=pmid)
            shutil.move('reach_output.json', prefix + '.json')
        else:
            rp = reach.process_json_file(prefix + '.json', citation=pmid)
        run_assembly(rp.statements, folder, pmcid)
