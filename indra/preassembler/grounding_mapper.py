import os
import csv
import sys
from copy import deepcopy
from indra.databases import uniprot_client
from itertools import groupby
from collections import Counter
import logging

logger = logging.getLogger('grounding_mapper')

class GroundingMapper(object):
    def __init__(self, gm):
        self.gm = gm

    def map_agents(self, stmts, do_rename=True):
        # Make a copy of the stmts
        mapped_stmts = []
        # Iterate over the statements
        for stmt in stmts:
            mapped_stmt = deepcopy(stmt)
            # Iterate over the agents
            skip_stmt = False
            for agent in mapped_stmt.agent_list():
                if agent is None or agent.db_refs.get('TEXT') is None:
                    continue
                agent_text = agent.db_refs.get('TEXT')
                # Look this string up in the grounding map
                # If not in the map, leave agent alone and continue
                if agent_text not in self.gm.keys():
                    pass
                # If it's in the map but it maps to None, then filter out
                # this statement by skipping it
                elif self.gm.get(agent_text) is None:
                    logger.debug("Skipping %s" % agent_text)
                    skip_stmt = True
                # If it has a value that's not None, map it and add it
                else:
                    # Otherwise, update the agent's db_refs field
                    agent.db_refs = self.gm.get(agent_text)
                    # Are we renaming right now?
                    if do_rename:
                        # If there's an INDRA ID, prefer that for the name
                        if agent.db_refs.get('INDRA'):
                            agent.name = agent.db_refs.get('INDRA')
                        # Take a HGNC name from Uniprot next
                        elif agent.db_refs.get('UP'):
                            # Try for the HGNC name
                            hgnc_name = uniprot_client.get_hgnc_name(
                                                      agent.db_refs.get('UP'))
                            if hgnc_name is not None:
                                agent.name = hgnc_name
                                agent.db_refs['HGNC'] = hgnc_name
                                continue
                            # Fall back on the Uniprot gene name
                            up_gene_name = uniprot_client.get_gene_name(
                                                       agent.db_refs.get('UP'))
                            if up_gene_name is not None:
                                agent.name = up_gene_name
                                continue
            # Check if we should skip the statement
            if not skip_stmt:
                mapped_stmts.append(mapped_stmt)
        return mapped_stmts

    def rename_agents(self, stmts):
        # Make a copy of the stmts
        mapped_stmts = deepcopy(stmts)
        # Iterate over the statements
        for stmt_ix, stmt in enumerate(mapped_stmts):
            # Iterate over the agents
            for agent in stmt.agent_list():
                if agent is None:
                    continue
                old_name = agent.name
                # If there's an INDRA ID, prefer that for the name
                if agent.db_refs.get('INDRA'):
                    agent.name = agent.db_refs.get('INDRA')
                # Take a HGNC name from Uniprot next
                elif agent.db_refs.get('UP'):
                    # Try for the HGNC name
                    hgnc_name = uniprot_client.get_hgnc_name(
                                                    agent.db_refs.get('UP'))
                    if hgnc_name is not None:
                        agent.name = hgnc_name
                        agent.db_refs['HGNC'] = hgnc_name
                        continue
                    # Fall back on the Uniprot gene name
                    up_gene_name = uniprot_client.get_gene_name(
                                                    agent.db_refs.get('UP'))
                    if up_gene_name is not None:
                        agent.name = up_gene_name
                        continue
                    # Take the text string
                    #if agent.db_refs.get('TEXT'):
                    #    agent.name = agent.db_refs.get('TEXT')
                    # If this fails, then we continue with no change
                # Fall back to the text string
                #elif agent.db_refs.get('TEXT'):
                #    agent.name = agent.db_refs.get('TEXT')
        return mapped_stmts

# TODO: handle the cases when there is more than one entry for the same
# key (e.g., ROS, ER)
def load_grounding_map(path):
    g_map = {}
    with open(path) as f:
        mapreader = csv.reader(f, delimiter=',', quotechar='"',
                               quoting=csv.QUOTE_MINIMAL,
                               lineterminator='\r\n')
        rows = [row for row in mapreader]

    for row in rows:
        key = row[0]
        db_refs = {'TEXT': key}
        keys = [entry for entry in row[1::2] if entry != '']
        values = [entry for entry in row[2::2] if entry != '']
        if len(keys) != len(values):
            logger.info('ERROR: Mismatched keys and values in row %s' %
                        str(row))
            continue
        else:
            db_refs.update(dict(zip(keys, values)))
            if len(db_refs.keys()) > 1:
                g_map[key] = db_refs
            else:
                g_map[key] = None
    return g_map


# Some useful functions for analyzing the grounding of sets of statements
# Put together all agent texts along with their grounding
def all_agents(stmts):
    agents = []
    for stmt in stmts:
        for agent in stmt.agent_list():
            if agent is not None:
                agents.append(agent)
    return agents


def agent_texts(agents):
    return [ag.db_refs.get('TEXT') for ag in agents]


def get_sentences_for_agent(text, stmts):
    sentences = []
    for stmt in stmts:
        for agent in stmt.agent_list():
            if agent is not None and agent.db_refs.get('TEXT') == text:
                sentences.append(stmt.evidence[0].text)
    return sentences


def agent_texts_with_grounding(stmts):
    allag = all_agents(stmts)
    for ag in allag:
        pfam_def = ag.db_refs.get('PFAM-DEF')
        if pfam_def is not None:
            ag.db_refs['PFAM-DEF'] = tuple(pfam_def)
    refs = [tuple(ag.db_refs.items()) for ag in allag]
    refs_counter = Counter(refs)
    refs_counter_dict = [(dict(entry[0]), entry[1])
                         for entry in refs_counter.items()]
    # First, sort by text
    refs_counter_dict_sorted = \
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
            if db_refs.keys() == ['TEXT']:
                db_ref_list.append((None, None, count))
            # Add any other db_refs (not TEXT)
            for db, id in db_refs.items():
                if db == 'TEXT':
                    continue
                else:
                    db_ref_list.append((db, id, count))
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
    ungrounded = [ag.db_refs['TEXT']
                  for s in stmts
                  for ag in s.agent_list()
                  if ag is not None and ag.db_refs.keys() == ['TEXT']]
    ungroundc = Counter(ungrounded)
    ungroundc = ungroundc.items()
    ungroundc.sort(key=lambda x: x[1], reverse=True)
    return ungroundc


def get_agents_with_name(name, stmts):
    return [ag for stmt in stmts for ag in stmt.agent_list()
               if ag is not None and ag.name == name]


def save_base_map(filename, grouped_by_text):
    rows = []
    for group in grouped_by_text:
        text_string = group[0]
        for db, id, count in group[1]:
            if db == 'UP':
                name = uniprot_client.get_mnemonic(id)
            else:
                name = ''
            row = [text_string.encode('utf8'), db, id, count, name]
            rows.append(row)

    with open(filename, 'w') as f:
        csvwriter = csv.writer(f, delimiter=',', quotechar='"',
                               quoting=csv.QUOTE_MINIMAL,
                               lineterminator='\r\n')
        csvwriter.writerows(rows)
        f.write('\r\n')

default_grounding_map_path = os.path.join(os.path.dirname(__file__),
                                  '../../bioentities/grounding_map.csv')
default_grounding_map = load_grounding_map(default_grounding_map_path)
gm = default_grounding_map


if __name__ == '__main__':
    import pickle

    if len(sys.argv) != 2:
        print "Usage: %s stmt_file" % sys.argv[0]
        sys.exit()
    statement_file = sys.argv[1]

    print "Opening statement file", statement_file
    with open(statement_file) as f:
        st = pickle.load(f)

    stmts = []
    for stmt_list in st.values():
        stmts += stmt_list

    twg = agent_texts_with_grounding(stmts)

    # Filter out those entries that are NOT already in the grounding map
    filtered_twg = [entry for entry in twg
                    if entry[0] not in default_grounding_map.keys()]

    save_base_map('%s_filtered_twg.csv' % statement_file, filtered_twg)

