import os
import csv
from collections import namedtuple
from copy import deepcopy
from indra.databases import uniprot_client
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
            # FIXME: Does not follow bound agents!!!
            # COMPLEXES
            # Map sites for complex
            if isinstance(stmt, Complex):
                invalid_sites = []
                stmt_copy.members = []
                for m in stmt.members:
                    (agent_invalid_sites, new_agent) = self.map_agent_sites(m)
                    stmt_copy.members.append(new_agent)
                    invalid_sites += agent_invalid_sites
                # If the list isn't empty, that means that there were incorrect
                # residues for this statement; add to mapped_statements list
                if invalid_sites:
                    mapped_stmt = \
                                MappedStatement(stmt, invalid_sites, stmt_copy)
                    mapped_statements.append(mapped_stmt)
                else:
                    valid_statements.append(stmt)
            # MODIFICATIONs
            # FIXME: Does not follow bound agents!!!
            elif isinstance(stmt, Modification):
                invalid_sites = []
                # Check substrate
                (agent_invalid_sites, new_sub) = self.map_agent_sites(stmt.sub)
                invalid_sites += agent_invalid_sites
                stmt_copy.sub = new_sub
                # Check enzyme
                (agent_invalid_sites, new_enz) = self.map_agent_sites(stmt.enz)
                invalid_sites += agent_invalid_sites
                stmt_copy.enz = new_enz

                # Check modification
                if stmt.residue is not None and stmt.position is not None:
                    assert isinstance(stmt.residue, basestring) and \
                           isinstance(stmt.position, basestring)
                    old_mod_list = [ModCondition(None, stmt.residue,
                                                 stmt.position)]
                    # Figure out if this site is invalid
                    stmt_invalid_sites = \
                            self.check_agent_mod(stmt_copy.sub, old_mod_list)
                    invalid_sites += stmt_invalid_sites
                    new_mod_list = \
                            update_mod_list(stmt_copy.sub.name, old_mod_list,
                                            stmt_invalid_sites)
                    stmt_copy.residue = new_mod_list[0].residue
                    stmt_copy.position = new_mod_list[0].position
                # Return valid/mapped site lists
                if invalid_sites:
                    mapped_stmt = \
                                MappedStatement(stmt, invalid_sites, stmt_copy)
                    mapped_statements.append(mapped_stmt)
                else:
                    valid_statements.append(stmt)
            """
            elif isinstance(stmt, ActivityModification):
                invalid_sites = []
                # Check agent
                (agent_invalid_sites, new_monomer) = \
                                self.map_agent_sites(stmt.monomer)
                invalid_sites += agent_invalid_sites
                stmt_copy.monomer = new_monomer
                # Check modification on sites
                # Filter lists
                if stmt.mod_pos and stmt.mod:
                    stmt_invalid_sites = \
                            self.check_agent_mod(stmt_copy.monomer,
                                                 filt_mod, filt_mod_pos)
                    invalid_sites += stmt_invalid_sites
                    (new_mod_list, new_modpos_list) = \
                            update_mod_list(stmt_copy.monomer.name,
                                            stmt_copy.mod,
                                            stmt_copy.mod_pos,
                                            stmt_invalid_sites)
                    stmt_copy.mod = new_mod_list
                    stmt_copy.modpos = new_modpos_list
                if invalid_sites:
                    mapped_stmt = \
                            MappedStatement(stmt, invalid_sites, stmt_copy)
                    mapped_statements.append(mapped_stmt)
                else:
                    valid_statements.append(stmt)
            """
        return (valid_statements, mapped_statements)

    """
    def check_sequence(self, stmt):
        #Check whether references to
        #residues and sequence positions are consistent with sequence
        #information in the UniProt database
        failures = []

        if isinstance(stmt, SelfModification):
            failures += check_agent_mod(stmt.sub)
            if stmt.mod_pos is not None:
                failures += \
                           check_agent_mod(stmt.enz, [stmt.mod], [stmt.mod_pos])
        elif isinstance(stmt, ActivityModification):
            failures += check_agent_mod(stmt.monomer)
            failures += check_agent_mod(stmt.monomer, stmt.mod, stmt.mod_pos)
        return failures
    """

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
        agent_entry = get_uniprot_entry(agent)
        # If the uniprot entry is not found, let it pass
        if not agent_entry:
            return invalid_sites # empty list
        # Look up all of the modifications in uniprot, and add them to the list
        # of invalid sites if they are missing
        for old_mod in mods:
            # If no site information for this residue, skip
            if old_mod.position is None or old_mod.residue is None:
                continue
            # Look up the residue/position in uniprot
            site_valid = uniprot_client.verify_location(agent_entry,
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


    def check_stmt_mod(self, agent, mods=None, mod_sites=None):
        """Look up the modification site in Uniprot and then the site map.
        """

        invalid_sites = {}
        # FIXME This should go outside the func?
        # If no UniProt ID is found, we don't report a failure
        up_id = agent.db_refs.get('UP')
        if up_id is None:
            return invalid_sites
        # FIXME
        # If the UniProt ID is a list then choose the first one.
        if not isinstance(up_id, basestring):
            up_id = up_id[0]
        agent_entry = uniprot_client.query_protein(up_id)

        # Figure out which modifications we are checking: the agent's
        # modifications or the modifications in a statement
        if mod_sites is not None:
            check_mods = mods
            check_mod_sites = mod_sites
        else:
            check_mods = agent.mods
            check_mod_sites = agent.mod_sites

        # Look up all of the modifications in uniprot, and add them to the list
        # of invalid sites if they are missing
        for m, mp in zip(check_mods, check_mod_sites):
            if mp is None:
                continue
            residue = site_abbrevs.get(m, None)
            if residue is None:
                continue
            ver = uniprot_client.verify_location(agent_entry, residue, mp)
            # If it's not found in Uniprot, then look it up in the site map
            if not ver:
                try:
                    mapped_site = \
                            self.site_map[(agent.name, str(residue), str(mp))]
                except KeyError:
                    mapped_site = None
                invalid_sites[(agent.name, str(residue), str(mp))] = mapped_site

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


def get_uniprot_entry(agent):
    # If no UniProt ID is found, we don't report a failure
    up_id = agent.db_refs.get('UP')
    if up_id is None:
        return None
    # If the UniProt ID is a list then choose the first one.
    if not isinstance(up_id, basestring):
        up_id = up_id[0]
    agent_entry = uniprot_client.query_protein(up_id)
    return agent_entry

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
                             'curated_site_map.txt')

default_site_map = load_site_map(default_site_map_path)

default_mapper = SiteMapper(default_site_map)
