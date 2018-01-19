


if __name__ == '__main__':
    import sys
    import boto3
    import botocore
    import time
    import pickle
    from indra.literature import pubmed_client
    from indra.tools.reading import submit_reading_pipeline as sub_aws
    from indra.tools import assemble_corpus as ac
    from indra.util import write_unicode_csv

    basename = sys.argv[1]
    # Get gene list
    with open('genes.txt', 'rt') as f:
        genes = [line.strip() for line in f.readlines()]

    # Assemble a list of PMIDs curated in Entrez gene
    pmids_for_genes = {}
    for gene_ix, gene in enumerate(genes):
        try:
            pmids = pubmed_client.get_ids_for_gene(gene)
        except ValueError:
            print("%s: Invalid gene name, skipping" % gene)
            continue
        print("%s: %d articles" % (gene, len(pmids)))
        pmids_for_genes[gene] = pmids
    pmids = set([pmid for pmid_list in pmids_for_genes.values()
                      for pmid in pmid_list])

    # Save the PMIDs to a file
    print("Saving PMIDs")
    with open('lab_meeting_pmids.txt', 'wt') as f:
        for pmid in pmids:
            f.write('%s\n' % pmid)

    #job_ids = sub_aws.submit_run_reach(basename, 'lab_meeting_pmids.txt',
    #                                   pmids_per_job=100000)
    #sub_aws.submit_combine(basename, job_ids)

    # Wait for the job to finish by checking S3
    client = boto3.client('s3')
    key = 'reading_results/%s/stmts/0_%d.pkl' % (basename, len(pmids))
    print("Looking for %s on S3" % key)
    while True:
        try:
            stmts_resp = client.get_object(Bucket='bigmech', Key=key)
            break
        except botocore.exceptions.ClientError as e:
            if e.response['Error']['Code'] =='NoSuchKey':
                print('Still processing...')
            # If there was some other kind of problem, re-raise the exception
        time.sleep(30)

    stmts_bytes = stmts_resp['Body'].read()
    stmts_by_paper = pickle.loads(stmts_bytes)
    stmts = [s for stmt_list in stmts_by_paper.values()
               for s in stmt_list]
    print("Grounding entities...")
    ground_stmts = ac.map_grounding(stmts)
    print("Detecting duplicate and overlapping statements...")
    stmts = ac.run_preassembly(ground_stmts)

    def get(agent_name, stmts):
        return [s for s in stmts if s.agent_list()[0] is not None
                and s.agent_list()[0].name == agent_name]

    lines = []
    for stmt in stmts:
        for ev in stmt.evidence:
            ag1 = ag2 = None
            if len(stmt.agent_list()) >= 1 and stmt.agent_list()[0]:
                ag1 = stmt.agent_list()[0].name
            if len(stmt.agent_list()) >= 2 and stmt.agent_list()[1]:
                ag2 = stmt.agent_list()[1].name
            line = [ag1, stmt.__class__.__name__, ag2,
                    str(stmt), ev.pmid, ev.text]
            lines.append(line)
    write_unicode_csv('%s_results.csv' % basename, lines)
