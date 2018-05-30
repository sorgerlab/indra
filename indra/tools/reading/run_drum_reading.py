import sys
import json
import time
import pickle
import logging
from indra.sources.trips import process_xml
from indra.sources.trips.drum_reader import DrumReader


logger = logging.getLogger('run_drum_reading')


def read_pmid_sentences(pmid_sentences, **drum_args):
    """Read sentences from a PMID-keyed dictonary and return all Statements

    Parameters
    ----------
    pmid_sentences : dict[str, list[str]]
        A dictonary where each key is a PMID pointing to a list of sentences
        to be read.

    **drum_args
        Keyword arguments passed directly to the DrumReader. Typical
        things to specify are `host` and `port`. If `run_drum` is specified
        as True, this process will internally run the DRUM reading system
        as a subprocess. Otherwise, DRUM is expected to be running
        independently.

    Returns
    -------
    all_statements : list[indra.statement.Statement]
        A list of INDRA Statements resulting from the reading
    """
    def _set_pmid(statements, pmid):
        for stmt in statements:
            for evidence in stmt.evidence:
                evidence.pmid = pmid

    # See if we need to start DRUM as a subprocess
    run_drum = drum_args.get('run_drum', False)
    drum_process = None
    all_statements = {}
    # Iterate over all the keys and sentences to read
    for pmid, sentences in pmid_sentences.items():
        logger.info('================================')
        logger.info('Processing %d sentences for %s' % (len(sentences), pmid))
        ts = time.time()
        # Make a DrumReader instance
        drum_args['name'] = 'DrumReader%s' % pmid
        dr = DrumReader(**drum_args)
        time.sleep(3)
        # If there is no DRUM process set yet, we get the one that was
        # just started by the DrumReader
        if run_drum and drum_process is None:
            drum_args.pop('run_drum', None)
            drum_process = dr.drum_system
            # By setting this, we ensuer that the reference to the
            # process is passed in to all future DrumReaders
            drum_args['drum_system'] = drum_process
        # Now read each sentence for this key
        for sentence in sentences:
            dr.read_text(sentence)
        # Start receiving results and exit when done
        try:
            dr.start()
        except SystemExit:
            pass
        statements = []
        # Process all the extractions into INDRA Statements
        for extraction in dr.extractions:
            # Sometimes we get nothing back
            if not extraction:
                continue
            tp = process_xml(extraction)
            statements += tp.statements
        # Set the PMIDs for the evidences of the Statements
        _set_pmid(statements, pmid)
        te = time.time()
        logger.info('Reading took %d seconds and produced %d Statements.' %
                    (te-ts, len(statements)))
        all_statements[pmid] = statements
    # If we were running a DRUM process, we should kill it
    if drum_process and dr.drum_system:
        dr._kill_drum()
    return all_statements


def read_text(text, **drum_args):
    """Read sentences from a PMID-keyed dictonary and return all Statements

    Parameters
    ----------
    text : str
        A block of text to run DRUM on

    **drum_args
        Keyword arguments passed directly to the DrumReader. Typical
        things to specify are 'host' and 'port'.

    Returns
    -------
    statements : list[indra.statement.Statement]
        A list of INDRA Statements resulting from the reading
    """
    return read_pmid_sentences({'PMID': text}, **drum_args)


def read_pmc(pmcid, **drum_args):
    # TODO: run DRUM in PMC reading mode here
    return

def save_results(statements, out_fname):
    with open(out_fname, 'wb') as fh:
        pickle.dump(statements, fh)


if __name__ == '__main__':
    file_name = sys.argv[0]
    host = sys.argv[1]
    port = sys.argv[2]
    with open(file_name, 'rt') as fh:
        content = json.load(fh)
    statements = read_pmid_sentences(content, host=host, port=port)
    save_results(statements, 'results.pkl')
