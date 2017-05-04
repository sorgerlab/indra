from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

if __name__ == '__main__':
    import sys
    import boto3
    import botocore.session
    from indra.literature import elsevier_client as ec

    bucket_name = 'bigmech'
    s3_client = boto3.client('s3')
    pmid_list_filename = sys.argv[1]
    job_name = sys.argv[2]
    pmid_list_key = 'pmid_lists/' + pmid_list_filename
    # Upload the pmid_list to Amazon S3
    s3_client.upload_file(pmid_list_filename, 'bigmech', pmid_list_key)

    # Get AWS credentials
    # http://stackoverflow.com/questions/36287720/boto3-get-credentials-dynamically
    session = botocore.session.get_session()
    access_key = session.get_credentials().access_key
    secret_key = session.get_credentials().secret_key

    # Get the Elsevier keys from the Elsevier client
    environment_vars = [
            {'name': ec.api_key_env_name,
             'value': ec.elsevier_keys.get('X-ELS-APIKey')},
            {'name': ec.inst_key_env_name,
             'value': ec.elsevier_keys.get('X-ELS-Insttoken')},
            {'name': 'AWS_ACCESS_KEY_ID',
             'value': access_key},
            {'name': 'AWS_SECRET_ACCESS_KEY',
             'value': secret_key}]
    # The command to run
    base_command_list = ['python', '-m',
                         'indra.tools.reading.run_reach_on_pmids_aws',
                         pmid_list_filename, '/tmp', '16']
    # Submit the reading job
    batch_client = boto3.client('batch')
    # Iterate over the list of PMIDs and submit the job in chunks
    batch_client.submit_job(jobName=job_name,
                            jobQueue='run_reach_queue',
                            jobDefinition='run_reach_jobdef',
                            containerOverrides={'environment': environment_vars,
                                                'command': command})

