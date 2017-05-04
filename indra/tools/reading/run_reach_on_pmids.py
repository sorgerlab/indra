from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import stat
import sys
import tempfile
import shutil
import subprocess
import glob
import json
import pickle
import logging
import functools
import multiprocessing as mp
from collections import Counter
from indra import reach
from indra.literature import pmc_client, s3_client, get_full_text, \
                             elsevier_client

# Logger
logger = logging.getLogger('runreach')

# Get a multiprocessing context
ctx = mp.get_context('spawn')

def join_parts(prefix):
    """Join different REACH output JSON files into a single JSON."""
    try:
        with open(prefix + '.uaz.entities.json', 'rt') as f:
            entities = json.load(f)
        with open(prefix + '.uaz.events.json', 'rt') as f:
            events = json.load(f)
        with open(prefix + '.uaz.sentences.json', 'rt') as f:
            sentences = json.load(f)
    except IOError as e:
        logger.error('Failed to open JSON files for %s; REACH error?' %
                      prefix)
        return None
    return {'events': events, 'entities': entities, 'sentences': sentences}


def upload_process_pmid(pmid_info, output_dir=None, reach_version=None):
    # The prefixes should be PMIDs
    pmid, source_text = pmid_info
    prefix_with_path = os.path.join(output_dir, pmid)
    full_json = join_parts(prefix_with_path)
    # Check that all parts of the JSON could be assembled
    if full_json is None:
        logger.error('REACH output missing JSON for %s' % pmid)
        return {pmid: []}
    # Upload the REACH output to S3
    s3_client.put_reach_output(full_json, pmid, reach_version, source_text)
    # Process the REACH output with INDRA
    return {pmid: process_reach_str(full_json, pmid)}


def upload_process_reach_files(output_dir, pmid_info_dict, reach_version,
                               num_cores):
    # At this point, we have a directory full of JSON files
    # Collect all the prefixes into a set, then iterate over the prefixes

    # Collect prefixes
    json_files = glob.glob(os.path.join(output_dir, '*.json'))
    json_prefixes = set([])
    for json_file in json_files:
        filename = os.path.basename(json_file)
        prefix = filename.split('.')[0]
        json_prefixes.add(prefix)
    # Make a list with PMID and source_text info
    pmid_info = [(json_prefix,
                  pmid_info_dict[json_prefix].get('content_source'))
                  for json_prefix in json_prefixes]
    # Create a multiprocessing pool
    logger.info('Creating a multiprocessing pool with %d cores' % num_cores)
    pool = ctx.Pool(num_cores)
    logger.info('Uploading and processing local REACH JSON files')
    upload_process_pmid_func = \
            functools.partial(upload_process_pmid, output_dir=output_dir,
                              reach_version=reach_version)
    res = pool.map(upload_process_pmid_func, pmid_info)

    """
    logger.info('Uploaded REACH JSON for %d files to S3 (%d failures)' %
        (num_uploaded, num_failures))
    failures_file = os.path.join(output_dir, 'failures.txt')
    with open(failures_file, 'wt') as f:
        for fail in failures:
            f.write('%s\n' % fail)
    """

# Version 1: If JSON is not available, get content and store;
#       assume force_read is False
# Version 1.5: If JSON is not available, get content and store;
#       check for force_read
# Version 2: If JSON is available, return JSON or process
# it and return statements (process it?)

def download_from_s3(pmid, input_dir=None, reach_version=None,
                     force_read=False, force_fulltext=False):
    if input_dir is None:
        raise ValueError('input_dir must be defined')

    # First define the text retrieval function
    def get_text():
        full_pmid = s3_client.check_pmid(pmid)
        # Look for the full text
        (content, content_type) = s3_client.get_upload_content(pmid,
                                        force_fulltext_lookup=force_fulltext)
        content_path = None
        # Write the contents to a file
        if content_type is None or content is None:
            # No content found on S3, skipping
            content_source = 'content_not_found'
        elif content_type == 'pmc_oa_xml':
            content_source = 'pmc_oa_xml'
            content_path = os.path.join(input_dir, '%s.nxml' % pmid)
        elif content_type == 'pmc_auth_xml':
            content_source = 'pmc_auth_xml'
            content_path = os.path.join(input_dir, '%s.nxml' % pmid)
        elif content_type == 'pmc_oa_txt':
            content_source = 'pmc_oa_txt'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'elsevier_xml':
            content = elsevier_client.extract_text(content)
            # Couldn't get text from Elsevier XML
            if content is None:
                content_source = 'elsevier_extract_text_failure'
            else:
                content_source = 'elsevier_xml'
                content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'txt':
            content_source = 'txt'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'abstract':
            content_source = 'abstract'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        # Unhandled content type, skipping
        else:
            content_source = 'unhandled_content_type_%s' % content_type
        # If we got content, write the content to a file with the appropriate
        # extension
        if content_path:
            with open(content_path, 'wb') as f:
                # The XML string is Unicode
                enc = content.encode('utf-8')
                f.write(enc)
        # Return dict of results for this PMID
        result = {pmid: {'content_source': content_source,
                         'content_path': content_path}}
        return result

    # If we're forcing a read regardless of whether there is cached REACH
    # output, then we download the text content
    if force_read or reach_version is None:
        return get_text()
    # If not, look for REACH JSON on S3
    (read_reach_version, read_source_text) = \
                            s3_client.get_reach_metadata(pmid)
    # Found it, same version, no need to get text
    if read_reach_version is not None and \
       read_reach_version == reach_version:
       result = {pmid: {'reach_version': read_reach_version,
                        'reach_source_text': read_source_text}}
    # Found it, different version, get the text
    else:
        result = get_text()
        result[pmid].update({'reach_version': read_reach_version,
                             'reach_source_text': read_source_text})
    return result


def process_reach_str(reach_json_str, pmid):
    if reach_json_str is None:
        raise ValueError('reach_json_str cannot be None')
    # Run the REACH processor on the JSON
    try:
        reach_proc = reach.process_json_str(reach_json_str, citation=pmid)
    # If there's a problem, skip it
    except Exception as e:
        print("Exception processing %s" % pmid)
        print(e)
        return []
    return reach_proc.statements


def process_reach_from_s3(pmid):
    reach_json_str = s3_client.get_reach_json_str(pmid)
    if reach_json_str is None:
        return []
    else:
        return {pmid: process_reach_str(reach_json_str, pmid)}


def run(pmid_list, tmp_dir, num_cores, start_index, end_index, force_read,
        force_fulltext, path_to_reach, reach_version, cleanup=False,
        verbose=True):
    if end_index > len(pmid_list):
        end_index = len(pmid_list)
    pmids_in_range = pmid_list[start_index:end_index]

    # Create the temp directories for input and output
    base_dir = tempfile.mkdtemp(prefix='read_%s_to_%s_' %
                                (start_index, end_index), dir=tmp_dir)
    # Make the temp directory writeable by REACH
    os.chmod(base_dir, stat.S_IRWXO | stat.S_IRWXU | stat.S_IRWXG)
    input_dir = os.path.join(base_dir, 'input')
    output_dir = os.path.join(base_dir, 'output')
    os.makedirs(input_dir)
    os.makedirs(output_dir)
    pmids_to_read = []

    # Get multiprocessing pool
    logger.info('Creating multiprocessing pool with %d cpus' % num_cores)
    pool = ctx.Pool(num_cores)
    logger.info('Getting content for PMIDs in parallel')
    download_from_s3_func = functools.partial(download_from_s3,
                                         input_dir=input_dir,
                                         reach_version=reach_version,
                                         force_read=force_read,
                                         force_fulltext=force_fulltext)
    res = pool.map(download_from_s3_func, pmids_in_range)
    # Close the pool to allow worker processes to exit before running REACH
    pool.close()
    # Combine the results into a single dict
    pmid_results = {pmid: results for pmid_dict in res
                                  for pmid, results in pmid_dict.items()}
    # Tabulate and log content results here
    pmids_read = {pmid: result for pmid, result in pmid_results.items()
                       if result.get('reach_version') == reach_version}
    pmids_unread = {pmid: pmid_results[pmid]
                    for pmid in
                    set(pmid_results.keys()).difference(set(pmids_read.keys()))}
    logger.info('%d / %d papers already read with REACH %s' %
                (len(pmids_read), len(pmid_results), reach_version))
    num_found = len([pmid for pmid in pmids_unread
                          if pmids_unread[pmid].get('content_path')])
    logger.info('Retrieved content for %d / %d papers to be read' %
                (num_found, len(pmids_unread)))
    # Tabulate sources and log in sorted order
    content_source_list = [pmid_dict.get('content_source')
                           for pmid_dict in pmids_unread.values()]
    content_source_counter = Counter(content_source_list)
    content_source_list = [(source, count)
                            for source, count in content_source_counter.items()]
    content_source_list.sort(key=lambda x: x[1], reverse=True)
    if content_source_list:
        logger.info('Content sources:')
        for source, count in content_source_list:
            logger.info('%s: %d' % (source, count))
    # Save text sources
    logger.info('Saving text sources...')
    text_source_file = os.path.join(base_dir, 'content_types.pkl')
    with open(text_source_file, 'wb') as f:
        pickle.dump(pmids_unread, f, protocol=2)
    import sys
    sys.exit()
    # Create the REACH configuration file
    conf_file_text = """
    #
    # Configuration file for reach
    #

    # this is the directory that stores the raw nxml, .csv, and/or .tsv files
    # this directory *must* exist
    papersDir = {input_dir}

    # this is where the output files containing the extracted mentions will be stored
    # if this directory doesn't exist it will be created
    outDir = {output_dir}

    # the output format for mentions: text, fries, indexcard, or assembly-csv (default is 'fries')
    outputTypes = ["fries"]

    # whether or not assembly should be run
    withAssembly = false

    # this is where the context files will be stored
    # if this directory doesn't exist it will be created
    contextDir = {output_dir}

    # this is where the brat standoff and text files are dumped
    bratDir = {output_dir}

    # verbose logging
    verbose = true

    # the encoding of input and output files
    encoding = "utf-8"


    # this is a list of sections that we should ignore
    ignoreSections = ["references", "materials", "materials|methods", "methods", "supplementary-material"]

    # context engine config
    contextEngine {{
        type = Policy4
        params = {{
            bound = 3
        }}
    }}

    logging {{
      # defines project-wide logging level
      loglevel = INFO
      logfile = {base_dir}/reach.log
    }}

    # restart configuration
    restart {{
      # restart allows batch jobs to skip over input files already successfully
      # processed
      useRestart = false
      # restart log is one filename per line list of input files already
      # successfully processed
      logfile = {base_dir}/restart.log
    }}

    # grounding configuration
    grounding: {{
      # List of AdHoc grounding files to insert, in order, into the grounding search sequence.
      # Each element of the list is a map of KB filename and optional meta info (not yet used):
      #   example: {{ kb: "adhoc.tsv", source: "NMZ at CMU" }}
      adHocFiles: [
        {{ kb: "NER-Grounding-Override.tsv.gz", source: "MITRE/NMZ/BG feedback overrides" }}
      ]

      # flag to turn off the influence of species on grounding
      overrideSpecies = true
    }}

    # number of simultaneous threads to use for parallelization
    threadLimit = {num_cores}

    # ReadPapers
    ReadPapers.papersDir = src/test/resources/inputs/nxml/
    ReadPapers.serializedPapers = mentions.ser
""".format(base_dir=base_dir, input_dir=input_dir, output_dir=output_dir,
               num_cores=num_cores)

    # Write the configuration file to the temp directory
    if num_found > 0:
        conf_file_path = os.path.join(base_dir, 'indra.conf')
        with open(conf_file_path, 'w') as f:
            f.write(conf_file_text)

        # Run REACH!
        args = ['java', '-jar', path_to_reach, conf_file_path]
        p = subprocess.Popen(args, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        if verbose:
            for line in iter(p.stdout.readline, b''):
                logger.info(line)
        (p_out, p_err) = p.communicate()
        if p.returncode:
            raise Exception(p_out.decode('utf-8') + '\n' +
                            p_err.decode('utf-8'))

        # Process JSON files from local file system, process to INDRA Statements
        # and upload to S3
        upload_process_reach_files(output_dir, pmids_unread, reach_version,
                                   num_cores)
        # Delete the tmp directory if desired
        if cleanup:
            shutil.rmtree(base_dir)

    # Create a new multiprocessing pool for processing the REACH JSON
    # files previously cached on S3
    logger.info('Creating multiprocessing pool with %d cpus' % num_cores)
    pool = ctx.Pool(num_cores)
    # Download and process the JSON files on S3
    logger.info('Processing REACH JSON from S3 in parallel')
    res = pool.map(process_reach_from_s3, pmids_read.keys())
    pool.close()


    # Save the list of PMIDs with no content found on S3/literature client
    #content_not_found_file = os.path.join(base_dir, 'content_not_found.txt')
    #with open(content_not_found_file, 'wt') as f:
    #    for c in content_not_found:
    #        f.write('%s\n' % c)



if __name__ == '__main__':
    # Set some variables
    cleanup = False
    verbose = True
    path_to_reach = '/pmc/reach/target/scala-2.11/reach-gordo-1.3.3-SNAPSHOT.jar'
    #path_to_reach = '/Users/johnbachman/Dropbox/1johndata/Knowledge File/Biology/Research/Big Mechanism/reach/target/scala-2.11/reach-gordo-1.3.3-SNAPSHOT.jar'
    reach_version = '1.3.3-b4a284'
    force_read = True
    force_fulltext = False

    # Check the arguments
    usage = "Usage: %s pmid_list tmp_dir num_cores start_index end_index " \
            "[force_fulltext]\n" % sys.argv[0]
    usage += "Alternative usage: %s upload_json output_dir content_types_file" % \
              sys.argv[0]
    if len(sys.argv) not in  (4, 6, 7):
        print(usage)
        sys.exit()
    if len(sys.argv) == 4 and sys.argv[1] != 'upload_json':
        print(usage)
        sys.exit()
    if len(sys.argv) == 7 and sys.argv[6] != 'force_fulltext':
        print(usage)
        sys.exit()
    elif len(sys.argv) == 7:
        force_fulltext = True

    # One type of operation: just upload previously read JSON files
    if len(sys.argv) == 4 and sys.argv[1] == 'upload_json':
        output_dir = sys.argv[2]
        text_sources_file = sys.argv[3]
        with open(text_sources_file, 'rb') as f:
            text_sources = pickle.load(f)
        upload_reach_json(output_dir, text_sources, reach_version)
        sys.exit()

    # =======================
    # Alternatively, run the whole process
    # Get the command line arguments
    (pmid_list_file, tmp_dir, num_cores, start_index, end_index) = sys.argv[1:6]
    start_index = int(start_index)
    end_index = int(end_index)
    num_cores = int(num_cores)

    # Load the list of PMIDs from the given file
    with open(pmid_list_file) as f:
        pmid_list = [line.strip('\n') for line in f.readlines()]

    # Do the reading
    run(pmid_list, tmp_dir, num_cores, start_index, end_index, force_read,
        force_fulltext, path_to_reach, reach_version, cleanup=cleanup,
        verbose=verbose)
