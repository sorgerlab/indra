from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import csv
import sys
import json
import pickle
import logging
from copy import deepcopy
from collections import Counter
from itertools import groupby, chain
from indra.statements import Agent
from indra.databases import uniprot_client, hgnc_client
from indra.util import read_unicode_csv, write_unicode_csv

logger = logging.getLogger('grounding_mapper')


class GroundingMapper(object):
    """Maps grounding of INDRA Agents based on a given grounding map.

    Attributes
    ----------
    gm : dict
        The grounding map, a dictionary mapping strings (entity names) to
        a dictionary of database identifiers.
    agent_map : Optional[dict]
        A dictionary mapping strings to grounded INDRA Agents with given state.
    """
    def __init__(self, gm, agent_map=None):
        self.gm = gm
        self.agent_map = agent_map if agent_map is not None else {}

    def update_agent_db_refs(self, agent, agent_text, do_rename=True):
        gene_name = None
        map_db_refs = deepcopy(self.gm.get(agent_text))
        up_id = map_db_refs.get('UP')
        hgnc_sym = map_db_refs.get('HGNC')
        if up_id and not hgnc_sym:
            gene_name = uniprot_client.get_gene_name(up_id, False)
            if gene_name:
                hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                if hgnc_id:
                    map_db_refs['HGNC'] = hgnc_id
        elif hgnc_sym and not up_id:
            # Override the HGNC symbol entry from the grounding
            # map with an HGNC ID
            hgnc_id = hgnc_client.get_hgnc_id(hgnc_sym)
            if hgnc_id:
                map_db_refs['HGNC'] = hgnc_id
                # Now get the Uniprot ID for the gene
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                if up_id:
                    map_db_refs['UP'] = up_id
            # If there's no HGNC ID for this symbol, raise an
            # Exception
            else:
                raise ValueError('No HGNC ID corresponding to gene '
                                 'symbol %s in grounding map.' %
                                 hgnc_sym)
        # If we have both, check the gene symbol ID against the
        # mapping from Uniprot
        elif up_id and hgnc_sym:
            # Get HGNC Symbol from Uniprot
            gene_name = uniprot_client.get_gene_name(up_id)
            if not gene_name:
                raise ValueError('No gene name found for Uniprot '
                                 'ID %s (expected %s)' %
                                 (up_id, hgnc_sym))
            # We got gene name, compare it to the HGNC name
            else:
                if gene_name != hgnc_sym:
                    raise ValueError('Gene name %s for Uniprot ID '
                                     '%s does not match HGNC '
                                     'symbol %s given in grounding '
                                     'map.' %
                                     (gene_name, up_id, hgnc_sym))
                else:
                    hgnc_id = hgnc_client.get_hgnc_id(hgnc_sym)
                    if not hgnc_id:
                        raise ValueError('No HGNC ID '
                                         'corresponding to gene '
                                         'symbol %s in grounding '
                                         'map.' % hgnc_sym)
                    map_db_refs['HGNC'] = hgnc_id
        # Assign the DB refs from the grounding map to the agent
        agent.db_refs = map_db_refs
        # Are we renaming right now?
        if do_rename:
            # If there's a FamPlex ID, prefer that for the name
            if agent.db_refs.get('FPLX'):
                agent.name = agent.db_refs.get('FPLX')
            # Get the HGNC symbol or gene name (retrieved above)
            elif hgnc_sym is not None:
                agent.name = hgnc_sym
            elif gene_name is not None:
                agent.name = gene_name
        return

    def map_agents_for_stmt(self, stmt, do_rename=True):
        mapped_stmt = deepcopy(stmt)
        # Iterate over the agents
        mapped_agent_list = mapped_stmt.agent_list()

        # Update agents directly participating in the statement
        agent_list = mapped_stmt.agent_list()
        for idx in range(len(agent_list)):
            agent = agent_list[idx]
            if agent is None or agent.db_refs.get('TEXT') is None:
                continue

            new_agent, maps_to_none = self.map_agent(agent, do_rename)

            if maps_to_none:
                # Skip the entire statement if the agent maps to None in the
                # grounding map
                return None

            # If the old agent had bound conditions, but the new agent does
            # not, copy the bound conditions over
            if new_agent is not None and len(new_agent.bound_conditions) == 0:
                new_agent.bound_conditions = agent.bound_conditions

            agent_list[idx] = new_agent

        mapped_stmt.set_agent_list(agent_list)

        # Update agents in the bound conditions
        for agent in agent_list:
            if agent is not None:
                for bc in agent.bound_conditions:
                    bc.agent, maps_to_none = self.map_agent(bc.agent, do_rename)
                    if maps_to_none:
                        # Skip the entire statement if the agent maps to None
                        # in the grounding map
                        return None

        return mapped_stmt

    def map_agent(self, agent, do_rename):
        """Grounds an agent; returns the new agent object (which might be
        a different object if we load a new agent state from javascript).

        Parameters
        ----------
        agent: indra.statements.Agent
            The agent to map
        do_rename: bool
            Whether to rename the agent text

        Returns
        -------
        grounded_agent: indra.statements.Agent
            The grounded agent
        maps_to_none: bool
            Whether the agent is in the grounding map and maps to None
        """

        agent_text = agent.db_refs.get('TEXT')
        mapped_to_agent_json = self.agent_map.get(agent_text)
        if mapped_to_agent_json:
            mapped_to_agent = \
                Agent._from_json(mapped_to_agent_json['agent'])
            return mapped_to_agent, False
        # Look this string up in the grounding map
        # If not in the map, leave agent alone and continue
        if agent_text in self.gm.keys():
            map_db_refs = self.gm[agent_text]
        else:
            return agent, False
        # If it's in the map but it maps to None, then filter out
        # this statement by skipping it
        if map_db_refs is None:
            # Increase counter if this statement has not already
            # been skipped via another agent
            logger.debug("Skipping %s" % agent_text)
            return None, True
        # If it has a value that's not None, map it and add it
        else:
            # Otherwise, update the agent's db_refs field
            self.update_agent_db_refs(agent, agent_text, do_rename)
        return agent, False

    def map_agents(self, stmts, do_rename=True):
        # Make a copy of the stmts
        mapped_stmts = []
        num_skipped = 0
        # Iterate over the statements
        for stmt in stmts:
            mapped_stmt = self.map_agents_for_stmt(stmt, do_rename)
            # Check if we should skip the statement
            if mapped_stmt is not None:
                mapped_stmts.append(mapped_stmt)
            else:
                num_skipped += 1
        logger.info('%s statements filtered out' % num_skipped)
        return mapped_stmts

    def rename_agents(self, stmts):
        # Make a copy of the stmts
        mapped_stmts = deepcopy(stmts)
        # Iterate over the statements
        for _, stmt in enumerate(mapped_stmts):
            # Iterate over the agents
            for agent in stmt.agent_list():
                if agent is None:
                    continue
                # If there's a FamPlex ID, prefer that for the name
                if agent.db_refs.get('FPLX'):
                    agent.name = agent.db_refs.get('FPLX')
                # Take a HGNC name from Uniprot next
                elif agent.db_refs.get('UP'):
                    # Try for the gene name
                    gene_name = uniprot_client.get_gene_name(
                                                    agent.db_refs.get('UP'),
                                                    web_fallback=False)
                    if gene_name:
                        agent.name = gene_name
                        hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                        if hgnc_id:
                            agent.db_refs['HGNC'] = hgnc_id
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
def load_grounding_map(grounding_map_path, ignore_path=None):
    g_map = {}
    map_rows = read_unicode_csv(grounding_map_path, delimiter=',',
                                quotechar='"',
                                quoting=csv.QUOTE_MINIMAL,
                                lineterminator='\r\n')
    if ignore_path and os.path.exists(ignore_path):
        ignore_rows = read_unicode_csv(ignore_path, delimiter=',',
                                       quotechar='"',
                                       quoting=csv.QUOTE_MINIMAL,
                                       lineterminator='\r\n')
    else:
        ignore_rows = []
    csv_rows = chain(map_rows, ignore_rows)
    for row in csv_rows:
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
            # Agents don't always have a TEXT db_refs entry (for instance
            # in the case of Statements from databases) so we check for this.
            if agent is not None and agent.db_refs.get('TEXT') is not None:
                agents.append(agent)
    return agents


def agent_texts(agents):
    return [ag.db_refs.get('TEXT') for ag in agents]


def get_sentences_for_agent(text, stmts, max_sentences=None):
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
    ungrounded = [ag.db_refs['TEXT']
                  for s in stmts
                  for ag in s.agent_list()
                  if ag is not None and list(ag.db_refs.keys()) == ['TEXT']]
    ungroundc = Counter(ungrounded)
    ungroundc = ungroundc.items()
    ungroundc = sorted(ungroundc, key=lambda x: x[1], reverse=True)
    return ungroundc


def get_agents_with_name(name, stmts):
    return [ag for stmt in stmts for ag in stmt.agent_list()
            if ag is not None and ag.name == name]


def save_base_map(filename, grouped_by_text):
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
    """Build map of entity texts to validated protein grounding.

    Looks at the grounding of the entity texts extracted from the statements
    and finds proteins where there is grounding to a human protein that maps to
    an HGNC name that is an exact match to the entity text. Returns a dict that
    can be used to update/expand the grounding map.
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


default_grounding_map_path = \
    os.path.join(os.path.dirname(__file__),
                 '../resources/famplex/grounding_map.csv')
default_ignore_path = \
    os.path.join(os.path.dirname(__file__),
                 '../resources/famplex/ignore.csv')
default_agent_grounding_path = \
    os.path.join(os.path.dirname(__file__),
                 '../resources/grounding_agents.json')
default_grounding_map = \
    load_grounding_map(default_grounding_map_path, default_ignore_path)


gm = default_grounding_map
with open(default_agent_grounding_path, 'r') as fh:
    default_agent_map = json.load(fh)


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print("Usage: %s stmt_file" % sys.argv[0])
        sys.exit()
    statement_file = sys.argv[1]

    logger.info("Opening statement file %s" % statement_file)
    with open(statement_file, 'rb') as f:
        st = pickle.load(f)

    stmts = []
    for stmt_list in st.values():
        stmts += stmt_list

    twg = agent_texts_with_grounding(stmts)

    save_base_map('%s_twg.csv' % statement_file, twg)

    # Filter out those entries that are NOT already in the grounding map
    filtered_twg = [entry for entry in twg
                    if entry[0] not in default_grounding_map.keys()]

    # For proteins that aren't explicitly grounded in the grounding map,
    # check for trivial corrections by building the protein map
    prot_map = protein_map_from_twg(twg)
    filtered_twg = [entry for entry in filtered_twg
                    if entry[0] not in prot_map.keys()]

    save_base_map('%s_unmapped_twg.csv' % statement_file, filtered_twg)

    # For each unmapped string, get sentences and write to file
    save_sentences(filtered_twg, stmts,
                   '%s_unmapped_sentences.csv' % statement_file)
