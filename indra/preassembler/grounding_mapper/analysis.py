import csv
import logging
from itertools import groupby
from collections import Counter
from indra.databases import uniprot_client
from indra.util import write_unicode_csv

logger = logging.getLogger(__name__)


# Some useful functions for analyzing the grounding of sets of statements
# Put together all agent texts along with their grounding
def all_agents(stmts):
    """Return a list of all of the agents from a list of statements.

    Only agents that are not None and have a TEXT entry are returned.

    Parameters
    ----------
    stmts : list of :py:class:`indra.statements.Statement`

    Returns
    -------
    agents : list of :py:class:`indra.statements.Agent`
        List of agents that appear in the input list of indra statements.
    """
    agents = []
    for stmt in stmts:
        for agent in stmt.agent_list():
            # Agents don't always have a TEXT db_refs entry (for instance
            # in the case of Statements from databases) so we check for this.
            if agent is not None and agent.db_refs.get('TEXT') is not None:
                agents.append(agent)
    return agents


def agent_texts(agents):
    """Return a list of all agent texts from a list of agents.

    None values are associated to agents without agent texts

    Parameters
    ----------
    agents : list of :py:class:`indra.statements.Agent`

    Returns
    -------
    list of str/None
        agent texts from input list of agents
    """
    return [ag.db_refs.get('TEXT') for ag in agents]


def get_sentences_for_agent(text, stmts, max_sentences=None):
    """Returns evidence sentences with a given agent text from a list of
    statements.

    Parameters
    ----------
    text : str
        An agent text

    stmts : list of :py:class:`indra.statements.Statement`
        INDRA Statements to search in for evidence statements.

    max_sentences : Optional[int/None]
        Cap on the number of evidence sentences to return. Default: None

    Returns
    -------
    sentences : list of str
        Evidence sentences from the list of statements containing
        the given agent text.
    """
    sentences = []
    for stmt in stmts:
        for agent in stmt.agent_list():
            if agent is not None and agent.db_refs.get('TEXT') == text:
                sentences.append((stmt.evidence[0].pmid,
                                  stmt.evidence[0].text))
                if max_sentences is not None and \
                        len(sentences) >= max_sentences:
                    return sentences
    return sentences


def agent_texts_with_grounding(stmts):
    """Return agent text groundings in a list of statements with their counts

    Parameters
    ----------
    stmts: list of :py:class:`indra.statements.Statement`

    Returns
    -------
    list of tuple
        List of tuples of the form
        (text: str, ((name_space: str, ID: str, count: int)...),
        total_count: int)

        Where the counts within the tuple of groundings give the number of
        times an agent with the given agent_text appears grounded with the
        particular name space and ID. The total_count gives the total number
        of times an agent with text appears in the list of statements.
    """
    allag = all_agents(stmts)
    # Convert PFAM-DEF lists into tuples so that they are hashable and can
    # be tabulated with a Counter
    for ag in allag:
        pfam_def = ag.db_refs.get('PFAM-DEF')
        if pfam_def is not None:
            ag.db_refs['PFAM-DEF'] = tuple(pfam_def)
    refs = [tuple(ag.db_refs.items()) for ag in allag]
    refs_counter = Counter(refs)
    refs_counter_dict = [(dict(entry[0]), entry[1])
                         for entry in refs_counter.items()]
    # First, sort by text so that we can do a groupby
    refs_counter_dict.sort(key=lambda x: x[0].get('TEXT'))

    # Then group by text
    grouped_by_text = []
    for k, g in groupby(refs_counter_dict, key=lambda x: x[0].get('TEXT')):
        # Total occurrences of this agent text
        total = 0
        entry = [k]
        db_ref_list = []
        for db_refs, count in g:
            # Check if TEXT is our only key, indicating no grounding
            if list(db_refs.keys()) == ['TEXT']:
                db_ref_list.append((None, None, count))
            # Add any other db_refs (not TEXT)
            for db, db_id in db_refs.items():
                if db == 'TEXT':
                    continue
                else:
                    db_ref_list.append((db, db_id, count))
            total += count
        # Sort the db_ref_list by the occurrences of each grounding
        entry.append(tuple(sorted(db_ref_list, key=lambda x: x[2],
                                  reverse=True)))
        # Now add the total frequency to the entry
        entry.append(total)
        # And add the entry to the overall list
        grouped_by_text.append(tuple(entry))
    # Sort the list by the total number of occurrences of each unique key
    grouped_by_text.sort(key=lambda x: x[2], reverse=True)
    return grouped_by_text


# List of all ungrounded entities by number of mentions
def ungrounded_texts(stmts):
    """Return a list of all ungrounded entities ordered by number of mentions

    Parameters
    ----------
    stmts : list of :py:class:`indra.statements.Statement`

    Returns
    -------
    ungroundc : list of tuple
       list of tuples of the form (text: str, count: int) sorted in descending
       order by count.
    """
    ungrounded = [ag.db_refs['TEXT']
                  for s in stmts
                  for ag in s.agent_list()
                  if ag is not None and list(ag.db_refs.keys()) == ['TEXT']]
    ungroundc = Counter(ungrounded)
    ungroundc = ungroundc.items()
    ungroundc = sorted(ungroundc, key=lambda x: x[1], reverse=True)
    return ungroundc


def get_agents_with_name(name, stmts):
    """Return all agents within a list of statements with a particular name."""
    return [ag for stmt in stmts for ag in stmt.agent_list()
            if ag is not None and ag.name == name]


def save_base_map(filename, grouped_by_text):
    """Dump a list of agents along with groundings and counts into a csv file

    Parameters
    ----------
    filename : str
        Filepath for output file
    grouped_by_text : list of tuple
        List of tuples of the form output by agent_texts_with_grounding
    """
    rows = []
    for group in grouped_by_text:
        text_string = group[0]
        for db, db_id, count in group[1]:
            if db == 'UP':
                name = uniprot_client.get_mnemonic(db_id)
            else:
                name = ''
            row = [text_string, db, db_id, count, name]
            rows.append(row)

    write_unicode_csv(filename, rows, delimiter=',', quotechar='"',
                      quoting=csv.QUOTE_MINIMAL, lineterminator='\r\n')


def protein_map_from_twg(twg):
    """Build  map of entity texts to validate protein grounding.

    Looks at the grounding of the entity texts extracted from the statements
    and finds proteins where there is grounding to a human protein that maps to
    an HGNC name that is an exact match to the entity text. Returns a dict that
    can be used to update/expand the grounding map.

    Parameters
    ----------
    twg : list of tuple
        list of tuples of the form output by agent_texts_with_grounding

    Returns
    -------
    protein_map : dict
        dict keyed on agent text with associated values
        {'TEXT': agent_text, 'UP': uniprot_id}. Entries are for agent texts
        where the grounding map was able to find human protein grounded to
        this agent_text in Uniprot.
    """

    protein_map = {}
    unmatched = 0
    matched = 0
    logger.info('Building grounding map for human proteins')
    for agent_text, grounding_list, _ in twg:
        # If 'UP' (Uniprot) not one of the grounding entries for this text,
        # then we skip it.
        if 'UP' not in [entry[0] for entry in grounding_list]:
            continue
        # Otherwise, collect all the Uniprot IDs for this protein.
        uniprot_ids = [entry[1] for entry in grounding_list
                       if entry[0] == 'UP']
        # For each Uniprot ID, look up the species
        for uniprot_id in uniprot_ids:
            # If it's not a human protein, skip it
            mnemonic = uniprot_client.get_mnemonic(uniprot_id)
            if mnemonic is None or not mnemonic.endswith('_HUMAN'):
                continue
            # Otherwise, look up the gene name in HGNC and match against the
            # agent text
            gene_name = uniprot_client.get_gene_name(uniprot_id)
            if gene_name is None:
                unmatched += 1
                continue
            if agent_text.upper() == gene_name.upper():
                matched += 1
                protein_map[agent_text] = {'TEXT': agent_text,
                                           'UP': uniprot_id}
            else:
                unmatched += 1
    logger.info('Exact matches for %d proteins' % matched)
    logger.info('No match (or no gene name) for %d proteins' % unmatched)
    return protein_map


def save_sentences(twg, stmts, filename, agent_limit=300):
    """Write evidence sentences for stmts with ungrounded agents to csv file.

    Parameters
    ----------
    twg: list of tuple
        list of tuples of ungrounded agent_texts with counts of the
        number of times they are mentioned in the list of statements.
        Should be sorted in descending order by the counts.
        This is of the form output by the function ungrounded texts.

    stmts: list of :py:class:`indra.statements.Statement`

    filename : str
        Path to output file

    agent_limit : Optional[int]
        Number of agents to include in output file. Takes the top agents
        by count.
    """
    sentences = []
    unmapped_texts = [t[0] for t in twg]
    counter = 0
    logger.info('Getting sentences for top %d unmapped agent texts.' %
                agent_limit)
    for text in unmapped_texts:
        agent_sentences = get_sentences_for_agent(text, stmts)
        sentences += map(lambda tup: (text,) + tup, agent_sentences)
        counter += 1
        if counter >= agent_limit:
            break
    # Write sentences to CSV file
    write_unicode_csv(filename, sentences, delimiter=',', quotechar='"',
                      quoting=csv.QUOTE_MINIMAL, lineterminator='\r\n')
