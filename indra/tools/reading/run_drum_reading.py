import sys
import json
import time
import pickle
from indra.sources.trips import process_xml
from indra.sources.trips.drum_reader import DrumReader


def read_pmid_sentences(pmid_sentences, **drum_args):
    """Read sentences from a PMID-keyed dictonary and return all Statements

    Parameters
    ----------
    pmid_sentences : dict[str, list[str]]
        A dictonary where each key is a PMID pointing to a list of sentences
        to be read.

    **drum_args
        Keyword arguments passed directly to the DrumReader. Typical
        things to specify are 'host' and 'port'.

    Returns
    -------
    all_statements : list[indra.statement.Statement]
        A list of INDRA Statements resulting from the reading
    """
    def _set_pmid(statements, pmid):
        for stmt in statements:
            for evidence in stmt.evidence:
                evidence.pmid = pmid

    all_statements = []
    for pmid, sentences in pmid_sentences.items():
        print('================================')
        print('Processing %d sentences for %s' % (len(sentences), pmid))
        ts = time.time()
        dr = DrumReader(**drum_args)
        for sentence in sentences:
            dr.read_text(sentence)
        try:
            dr.start()
        except SystemExit:
            pass
        statements = []
        for extraction in dr.extractions:
            tp = process_xml(extraction)
            statements += tp.statements
        _set_pmid(statements, pmid)
        te = time.time()
        print('Reading took %d seconds and produced %d Statements.' %
              (te-ts, len(statements)))
        all_statements += statements
    return all_statements


def read_text(text, **drum_args):
    """Read sentences from a PMID-keyed dictonary and return all Statements

    Parameters
    ----------
    pmid_sentences : dict[str, list[str]]
        A dictonary where each key is a PMID pointing to a list of sentences
        to be read.

    **drum_args
        Keyword arguments passed directly to the DrumReader. Typical
        things to specify are 'host' and 'port'.

    Returns
    -------
    all_statements : list[indra.statement.Statement]
        A list of INDRA Statements resulting from the reading
    """



def read_pmc(pmcid, **drum_args):
    


def save_results(statements, out_fname):
    with open(out_fname, 'wb') as fh:
        pickle.dump(statements, fh)


if __name__ == '__main__':
    host = sys.argv[1]
    file_name = sys.argv[2]
    with open(file_name, 'rt') as fh:
        content = json.load(fh)
    statements = read_content(content, host=host, port=port)
    save_results(statements, 'results.pkl')
