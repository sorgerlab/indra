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
from indra.literature import pmc_client, s3_client, get_full_text, \
                             elsevier_client

# Logger
logger = logging.getLogger('runreach')

def upload_reach_json(output_dir, text_sources):
    # At this point, we have a directory full of JSON files
    # Collect all the prefixes into a set, then iterate over the prefixes
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

    # Collect prefixes
    json_files = glob.glob(os.path.join(output_dir, '*.json'))
    json_prefixes = set([])
    for json_file in json_files:
        filename = os.path.basename(json_file)
        prefix = filename.split('.')[0]
        json_prefixes.add(prefix)
    # Now iterate over the collected prefixes, combine the JSON, and send to S3
    num_uploaded = 0
    num_failures = 0
    failures = []
    # The prefixes should be PMIDs
    for json_ix, json_prefix in enumerate(json_prefixes):
        prefix_with_path = os.path.join(output_dir, json_prefix)
        full_json = join_parts(prefix_with_path)
        if full_json is None:
            num_failures += 1
            failures.append(json_prefix)
        else:
            # Look up the paper source type
            source_text = text_sources.get(json_prefix)
            logger.info('%s (%d of %d): source %s' %
                      (json_prefix, json_ix + 1, len(json_prefixes), source_text))
            s3_client.put_reach_output(full_json, json_prefix, reach_version,
                                       source_text)
            num_uploaded += 1
    logger.info('Uploaded REACH JSON for %d files to S3 (%d failures)' %
        (num_uploaded, num_failures))
    failures_file = os.path.join(output_dir, 'failures.txt')
    with open(failures_file, 'wt') as f:
        for fail in failures:
            f.write('%s\n' % fail)


if __name__ == '__main__':

    cleanup = False
    verbose = True
    path_to_reach = '/pmc/reach/target/scala-2.11/reach-gordo-1.3.3-SNAPSHOT.jar'
    #path_to_reach = '/Users/johnbachman/Dropbox/1johndata/Knowledge File/Biology/Research/Big Mechanism/reach/target/scala-2.11/reach-gordo-1.3.3-SNAPSHOT.jar'
    reach_version = '1.3.3-38d2c4330'
    force_read = False
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
        upload_reach_json(output_dir, text_sources)
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

    # If we're re-reading no matter what, we don't have to check for existing
    # REACH output
    if force_read:
        pmids_to_read = pmids_in_range
    # Otherwise, check if we've read the PMIDs already
    else:
        for pmid in pmids_in_range:
            pmid = s3_client.check_pmid(pmid)
            (read_reach_version, read_source_text) = \
                                    s3_client.get_reach_metadata(pmid)
            # Found it, same version
            if read_reach_version is not None and \
               read_reach_version == reach_version:
                logger.info('%s: found same version (%s), skipping' %
                            (pmid, read_reach_version))
            # Found it, different version
            else:
                logger.info('%s: found %s, current %s; will re-read' %
                            (pmid, read_reach_version, reach_version))
                pmids_to_read.append(pmid)

    if not pmids_to_read:
        logger.info('No pmids to read!')
        sys.exit(0)
    else:
        logger.info('Preparing to run REACH on %d PMIDs' % len(pmids_to_read))

    # Now iterate over the pmids to read and download from S3 to the input
    # directory
    num_pmc_oa_xml = 0
    num_pmc_auth_xml = 0
    num_txt = 0
    num_elsevier_xml = 0
    num_abstract = 0
    num_not_found = 0
    num_elsevier_xml_fail = 0
    # Keep a map of the content type we've downloaded for each PMID
    text_sources = {}
    content_not_found = []
    for pmid in pmids_to_read:
        full_pmid = s3_client.check_pmid(pmid)
        # Look for the full text
        (content, content_type) = s3_client.get_upload_content(pmid, force_fulltext_lookup=force_fulltext)
        # If we don't find the XML on S3, look for it using the PMC client
        #if xml:
        #    num_found_s3 += 1
        #else:
        #    logger.info('No content for %s from S3' % pmid)
        #    (content, content_type) = get_full_text(pmid, 'pmid')
        #    if content_type == 'nxml':
        #        logger.info('Found nxml for %s from PMC web service' % pmid)
        #        xml = content
        #        num_found_not_s3 += 1
        #        # Upload the xml to S3 for next time
        #        logger.info('Uploading full text for %s to S3' % pmid)
        #        s3_client.put_full_text(pmid, xml, full_text_type='pmc_oa_xml')
        #    #elif content_type == 'abstract':
        #    #    logger.info('Found abstract for %s' % pmid)
        #    #    s3_client.put_abstract(pmid, content)
        #    else:
        #        logger.info('No full text found for %s' % pmid)
        # Write the contents to a file
        if content_type is None or content is None:
            num_not_found += 1
            content_not_found.append(pmid)
            logger.info('No content found on S3 for %s, skipping' % pmid)
            continue
        elif content_type == 'pmc_oa_xml':
            num_pmc_oa_xml += 1
            text_sources[full_pmid] = 'pmc_oa_xml'
            content_path = os.path.join(input_dir, '%s.nxml' % pmid)
        elif content_type == 'pmc_auth_xml':
            num_pmc_auth_xml += 1
            text_sources[full_pmid] = 'pmc_auth_xml'
            content_path = os.path.join(input_dir, '%s.nxml' % pmid)
        elif content_type == 'pmc_oa_txt':
            num_txt += 1
            text_sources[full_pmid] = 'pmc_oa_txt'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'elsevier_xml':
            content = elsevier_client.extract_text(content)
            if content is None:
                logger.info("%s: Couldn't get text from Elsevier XML" % pmid)
                num_elsevier_xml_fail += 1
                continue
            num_elsevier_xml += 1
            text_sources[full_pmid] = 'elsevier_xml'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'txt':
            num_txt += 1
            text_sources[full_pmid] = 'txt'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'abstract':
            num_abstract += 1
            text_sources[full_pmid] = 'abstract'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        else:
            num_not_found += 1
            logger.info('Unhandled content type %s for %s, skipping' %
                        (content_type, pmid))
            continue
        # Write the content to a file with the appropriate extension
        with open(content_path, 'wb') as f:
            # The XML string is Unicode
            enc = content.encode('utf-8')
            f.write(enc)
    logger.info('Saving text sources...')
    text_source_file = os.path.join(base_dir, 'content_types.pkl')
    with open(text_source_file, 'wb') as f:
        pickle.dump(text_sources, f, protocol=2)
    logger.info('Found content PMIDs:')
    logger.info('%d pmc_oa_xml' % num_pmc_oa_xml)
    logger.info('%d pmc_auth_xml' % num_pmc_auth_xml)
    logger.info('%d elsevier_xml' % num_elsevier_xml)
    logger.info('%d elsevier_xml with no full text' % num_elsevier_xml_fail)
    logger.info('%d txt (incl. some Elsevier)' % num_txt)
    logger.info('%d abstract' % num_abstract)
    logger.info('%d no content' % num_not_found)

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
    conf_file_path = os.path.join(base_dir, 'indra.conf')
    with open(conf_file_path, 'w') as f:
        f.write(conf_file_text)

    # Run REACH!
    args = ['java', '-jar', path_to_reach, conf_file_path]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if verbose:
        for line in iter(p.stdout.readline, b''):
            logger.info(line)
    (p_out, p_err) = p.communicate()
    if p.returncode:
        raise Exception(p_out.decode('utf-8') + '\n' + p_err.decode('utf-8'))

    # Save the list of PMIDs with no content found on S3/literature client
    content_not_found_file = os.path.join(base_dir, 'content_not_found.txt')
    with open(content_not_found_file, 'wt') as f:
        for c in content_not_found:
            f.write('%s\n' % c)

    # Upload!
    upload_reach_json(output_dir, text_sources)

    if cleanup:
        shutil.rmtree(base_dir)


