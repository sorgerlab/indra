import os
import csv
from copy import deepcopy
from indra.databases import uniprot_client
from itertools import groupby
from collections import Counter
import logging

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
                    print "Skipping %s" % agent_text
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
                if old_name != agent.name:
                    print "Map %d of %d: %s (%s) --> %s" % \
                                (stmt_ix+1, len(stmts), old_name,
                                 agent.db_refs.get('TEXT'), agent.name)
        return mapped_stmts

# TODO: handle the cases when there is more than one entry for the same
# key (e.g., ROS, ER)
def load_grounding_map(path):
    g_map = {}
    with open(path) as f:
        mapreader = csv.reader(f, delimiter='\t')
        for row in mapreader:
            key = row[0]
            db_refs = {'TEXT': key}
            for pair_ix in range(0, 2):
                col_ix = (pair_ix * 2) + 1
                db = row[col_ix]
                db_id = row[col_ix + 1]
                if db == '' or db == 'None' or db_id == '' or db_id == 'None':
                    continue
                else:
                    db_refs[db] = db_id
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
    refs_counter_dict_sorted = \
            refs_counter_dict.sort(key=lambda x: x[0].get('TEXT'))

    grouped_by_text = []
    for k, g in groupby(refs_counter_dict, key=lambda x: x[0].get('TEXT')):
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
        entry.append(tuple(sorted(db_ref_list, key=lambda x: x[2],
                     reverse=True)))
        entry.append(total)
        grouped_by_text.append(tuple(entry))

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
    with open(filename, 'w') as f:
        for group in grouped_by_text:
            text_string = group[0]
            for db, id, count in group[1]:
                if db == 'UP':
                    name = uniprot_client.get_mnemonic(id)
                else:
                    name = ''
                line = '%s\t%s\t%s\t%s\t%s\n' % \
                        (text_string, db, id, count, name)
                f.write(line.encode('utf8'))


default_grounding_map_path = os.path.join(os.path.dirname(__file__),
                                  '../resources/grounding_map_curated.txt')
default_grounding_map = load_grounding_map(default_grounding_map_path)
gm = default_grounding_map


if __name__ == '__main__':
    import pickle

    with open('trips_b4_stmts.pkl') as f:
        st = pickle.load(f)

    stmts = []
    for stmt_list in st.values():
        stmts += stmt_list

    twg = agent_texts_with_grounding(stmts)

    # Filter out those entries that are NOT already in the grounding map
    filtered_twg = [entry for entry in twg
                    if entry[0] not in default_grounding_map.keys()]

    #save_base_map('trips_batch4_base_map.txt', filtered_twg)

