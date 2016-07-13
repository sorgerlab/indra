import sys
import shutil
from indra import reach
from indra.literature import pmc_client, get_full_text, id_lookup
from assembly_eval import have_file, run_assembly

if __name__ == '__main__':
    folder = 'reach'
    pmc_ids = [s.strip() for s in open('pmcids.txt', 'rt').readlines()]
    pmids = [id_lookup(pmcid)['pmid'] for pmcid in pmc_ids]
    # Set to True only if reading should be ran again
    rerun = False

    # Download the papers if they are not available yet
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

    # Read each paper if it hasn't been read yet.
    # Otherwise use the existing json extractions.
    for pmcid, pmid in zip(pmc_ids, pmids):
        prefix = folder + '/' + pmcid
        print 'Processing %s...' % pmcid
        # If REACH already processed it then don't run it again
        if rerun or not have_file(prefix + '.json'):
            if have_file(prefix + '.txt'):
                txt = open(prefix + '.txt').read().decode('utf-8')
                rp = reach.process_text(txt, citation=pmid, offline=True)
            elif have_file(prefix + '.nxml'):
                print 'FOUND'
                rp = reach.process_nxml_file(prefix + '.nxml', citation=pmid,
                                             offline=True)
            if rp is None:
                print 'Reading with offline REACH failed.'
                continue
            shutil.move('reach_output.json', prefix + '.json')
        else:
            rp = reach.process_json_file(prefix + '.json', citation=pmid)
        run_assembly(rp.statements, folder, pmcid)
