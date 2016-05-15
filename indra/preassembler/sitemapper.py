import os
import csv
import warnings
from collections import namedtuple
from copy import deepcopy
from indra.databases import uniprot_client, hgnc_client
from indra.statements import *


MappedStatement = namedtuple('MappedStatement',
                             ['original_stmt', 'mapped_mods', 'mapped_stmt'])


class SiteMapper(object):
    """
    Parameters
    ----------
    site_map : dict mapping tuples to tuples as returned by `load_site_map`.
        A dict mapping tuples of the form `(gene, orig_res, orig_pos)` to a
        tuple of the form `(correct_res, correct_pos, comment)`, where `gene`
        is the string name of the gene (canonicalized to HGNC); `orig_res` and
        `orig_pos` are the residue and position to be mapped; `correct_res` and
        `correct_pos` are the corrected residue and position, and `comment` is
        a string describing the reason for the mapping (species error, isoform
        error, wrong residue name, etc.).
    """

    def __init__(self, site_map):
        self.site_map = site_map

    def map_sites(self, stmts, save_fname=None):
        """Iterates over a list of statements and runs checks on them.  Then it
        returns a tuple of lists, with the first element containing statements
        that passed all checks, and the second the statements that failed the
        tests

        If there is nothing amiss with the statement (modifications on any
        of the agents, modifications made in the statement, etc, then the
        statement goes into the valid_stmts list.

        If there is a problem with the statement, the offending modifications
        are looked up in the site map. If corresponding entries are found,
        a MappedStatement instance is created.
        """

        valid_statements = []
        mapped_statements = []

        for stmt in stmts:
            stmt_copy = deepcopy(stmt)
            # For all statements, replace agents with invalid modifications
            invalid_sites = []
            new_agent_list = []
            for agent in stmt.agent_list():
                (agent_invalid_sites, new_agent) = self.map_agent_sites(agent)
                new_agent_list.append(new_agent)
                invalid_sites += agent_invalid_sites
            if invalid_sites:
                stmt_copy.set_agent_list(new_agent_list)

            # --- Special handling for these statements ---
            # For modifications, fix residue and position
            if (isinstance(stmt, Modification) or \
                isinstance(stmt, SelfModification)) and \
                 stmt.residue is not None and stmt.position is not None:
                # Make sure we didn't end up with lists by accident
                assert isinstance(stmt.residue, basestring) and \
                       isinstance(stmt.position, basestring)
                # Get the right agent depending on whether this is a
                # Modification or SelfModification statement
                agent_to_check = (stmt_copy.sub
                                  if isinstance(stmt, Modification)
                                  else stmt_copy.enz)
                # Check the modification on the appropriate agent
                old_mod_list = [ModCondition(None, stmt.residue,
                                             stmt.position)]
                # Figure out if this site is invalid
                stmt_invalid_sites = \
                        self.check_agent_mod(agent_to_check, old_mod_list)
                # Add to our list of invalid sites
                invalid_sites += stmt_invalid_sites
                # Get the updated list of ModCondition objects
                new_mod_list = \
                        update_mod_list(agent_to_check.name, old_mod_list,
                                        stmt_invalid_sites)
                # Update the statement with the correct site
                stmt_copy.residue = new_mod_list[0].residue
                stmt_copy.position = new_mod_list[0].position

            # If the invalid_sites list isneck_agent_mod(self, agent, mods)t empty, that means that there were
            # incorrect residues for this statement; add the statement to
            # the mapped_statements list
            if invalid_sites:
                mapped_stmt = \
                            MappedStatement(stmt, invalid_sites, stmt_copy)
                mapped_statements.append(mapped_stmt)
            else:
                valid_statements.append(stmt)

        return (valid_statements, mapped_statements)

    def map_agent_sites(self, agent):
        """Look up the modification site in Uniprot and then the site map.
        """
        new_agent = deepcopy(agent)
        if not agent.mods:
            return ([], new_agent)
        invalid_sites = self.check_agent_mod(agent, agent.mods)
        if not invalid_sites:
            return ([], new_agent)
        new_mod_list = update_mod_list(agent.name, agent.mods, invalid_sites)
        # Finally, update the agent, and return along with invalid site info
        new_agent.mods = new_mod_list
        return (invalid_sites, new_agent)

    def check_agent_mod(self, agent, mods):
        """Look up the modification site in Uniprot and then the site map.
        Return a list of invalid sites, where each entry in the list has
        two elements: ((gene_name, residue, position), mapped_site).
        If the invalid position was not found in the site map, mapped_site is
        None; otherwise it is a tuple consisting of (residue, position,
        comment).
        """

        invalid_sites = []
        up_id = get_uniprot_id(agent)
        # If the uniprot entry is not found, let it pass
        if not up_id:
            return invalid_sites # empty list
        # Look up all of the modifications in uniprot, and add them to the list
        # of invalid sites if they are missing
        for old_mod in mods:
            # If no site information for this residue, skip
            if old_mod.position is None or old_mod.residue is None:
                continue
            # Look up the residue/position in uniprot
            site_valid = uniprot_client.verify_location(up_id,
                                                        old_mod.residue,
                                                        old_mod.position)
            # If it's not found in Uniprot, then look it up in the site map
            if not site_valid:
                site_key = (agent.name, old_mod.residue, old_mod.position)
                mapped_site = self.site_map.get(site_key, None)
                # We found an entry in the site map!
                if mapped_site is not None:
                    invalid_sites.append((site_key, mapped_site))
                # No entry in the site map--set site info to None
                else:
                    invalid_sites.append((site_key, None))
        return invalid_sites


def update_mod_list(agent_name, mods, invalid_sites):
    new_mod_list = []
    # Get the list of invalid/mapped sites for the agent
    invalid_site_keys = [site[0] for site in invalid_sites]
    for old_mod in mods:
        old_mod_key = (agent_name, old_mod.residue, old_mod.position)
        # If the original modification was found to be invalid, create a newly
        # updated modification
        if old_mod_key in invalid_site_keys:
            mapped_site = \
                    invalid_sites[invalid_site_keys.index(old_mod_key)][1]
            # No entry in the map: pass the incorrect site through
            if mapped_site is None:
                new_mod_list.append(old_mod)
            # Entry in the map
            else:
                # Do we have actual site information?
                new_res = mapped_site[0]
                new_pos = mapped_site[1]
                if new_res is not None and new_pos is not None:
                    new_mod_list.append(
                            ModCondition(old_mod.mod_type, new_res, new_pos,
                                         old_mod.is_modified))
                # Mapped, but no site info--pass through unchanged
                else:
                    new_mod_list.append(old_mod)
        # The modification is not in the invalid site list, so it's considered
        # valid
        else:
            new_mod_list.append(old_mod)
    return new_mod_list


def get_uniprot_id(agent):
    up_id = agent.db_refs.get('UP')
    hgnc_id = agent.db_refs.get('HGNC')
    if up_id is None:
        if hgnc_id is None:
            # If both UniProt and HGNC refs are missing we can't
            # sequence check and so don't report a failure.
            return None
        # Try to get UniProt ID from HGNC
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        # If this fails, again, we can't sequence check
        if up_id is None:
            return None

    # If the UniProt ID is a list then choose the first one.
    if not isinstance(up_id, basestring) and \
       isinstance(up_id[0], basestring):
        up_id = up_id[0]
    return up_id


def load_site_map(path):
    """Load the modification site map from a file.

    Parameters
    ----------
    path : string
        Path to the tab-separated site map file.

    Returns
    -------
    A dict mapping tuples of the form `(gene, orig_res, orig_pos)` to
    a tuple of the form `(correct_res, correct_pos, comment)`, where
    `gene` is the string name of the gene (canonicalized to HGNC); `orig_res`
    and `orig_pos` are the residue and position to be mapped; `correct_res`
    and `correct_pos` are the corrected residue and position, and `comment` is
    a string describing the reason for the mapping (species error, isoform
    error, wrong residue name, etc.).
    """

    site_map = {}
    with open(path) as f:
        mapreader = csv.reader(f, delimiter='\t')
        # Skip the header line
        mapreader.next()
        for row in mapreader:
            # Don't allow empty entries in the key section
            if not (row[0] and row[1] and row[2]):
                raise Exception("Entries in the key (gene, residue, position) "
                                "may not be empty.")
            correct_res = row[3].strip() if row[3] else None
            correct_pos = row[4].strip() if row[4] else None
            comment = row[5].strip() if row[5] else None
            site_map[(row[0].strip(), row[1].strip(), row[2].strip())] = \
                                    (correct_res, correct_pos, comment)
    return site_map


default_site_map_path = os.path.join(os.path.dirname(__file__),
                             '../resources/curated_site_map.txt')

default_site_map = load_site_map(default_site_map_path)

default_mapper = SiteMapper(default_site_map)
