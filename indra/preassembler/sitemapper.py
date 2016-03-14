import os
import csv
from collections import namedtuple
from copy import deepcopy
from indra.pysb_assembler import abbrevs as site_abbrevs
from indra.databases import uniprot_client
from indra.statements import *

MappedStatement = namedtuple('MappedStatement',
                             ['original_stmt', 'mapped_mods', 'mapped_stmt'])
UnmappedStatement = namedtuple('UnmappedStatement',
                               ['original_stmt', 'unmapped_mods',
                                'is_curated', 'comment'])
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

        If corresponding entres are not found for all residues, an
        UnmappedStatement entry is created.
        """

        pass_stmts = []
        fail_stmts = []
        failures = []

        for stmt in stmts:
            result = check_sequence(stmt)
            #failures += check_sequence(stmt)
            #if failures:
            #    fail_stmts.append(stmt)
            #else:
            #    pass_stmts.append(stmt)

    def check_sequence(self, stmt):
        """Check whether references to
        residues and sequence positions are consistent with sequence
        information in the UniProt database"""
        failures = []

        # Map sites for complex
        stmt_copy = deepcopy(stmt)
        if isinstance(stmt, Complex):
            for m in stmt.members:
                agent_mod_result = check_agent_mod(m)
                # Check the list of invalid sites--if it is not the empty
                # dict, then we replace this agent in the new statement with
                # the mapped agent

        elif isinstance(stmt, Modification):
            failures += check_agent_mod(stmt.sub)
            failures += check_agent_mod(stmt.enz)
            if stmt.mod_pos is not None:
                failures += check_agent_mod(stmt.sub, [stmt.mod], [stmt.mod_pos])
        elif isinstance(stmt, SelfModification):
            failures += check_agent_mod(stmt.sub)
            if stmt.mod_pos is not None:
                failures += check_agent_mod(stmt.enz, [stmt.mod], [stmt.mod_pos])
        elif isinstance(stmt, ActivityModification):
            failures += check_agent_mod(stmt.monomer)
            failures += check_agent_mod(stmt.monomer, stmt.mod, stmt.mod_pos)
        return failures

    def check_agent_mod(self, agent):
        """Look up the modification site in Uniprot and then the site map.
        """
        invalid_sites = {}
        new_agent = deepcopy(agent)
        agent_entry = get_uniprot_entry(agent)
        # If the uniprot entry is not found, let it pass
        if not agent_entry:
            return ({}, agent)

        # Look up all of the modifications in uniprot, and add them to the list
        # of invalid sites if they are missing
        new_mod_list = []
        new_modpos_list = []
        for old_mod, old_modpos in zip(agent.mods, agent.mod_sites):
            # If the modification doesn't have a site, add it to the list
            # without doing any lookup
            if old_mod is None:
                new_mod_list.append(old_mod)
                new_modpos_list.append(old_modpos)
                continue
            # Get the amino acid abbreviation (e.g., 'S', 'T', 'Y')
            residue = site_abbrevs.get(old_mod, None)
            # If no site for this type of modification, add it to agent list
            # and continue
            if residue is None:
                new_mod_list.append(old_mod)
                new_modpos_list.append(old_modpos)
                continue
            # Look up the residue/position in uniprot
            ver = uniprot_client.verify_location(agent_entry, residue,
                                                 old_modpos)
            # If it's not found in Uniprot, then look it up in the site map
            if not ver:
                # Look the site up in the sitemap
                mapped_site = self.site_map.get(
                             (agent.name, str(residue), str(old_modpos)), None)
                # We found an entry in the site map
                if mapped_site is not None:
                    new_res = mapped_site[0]
                    new_pos = mapped_site[1]
                    # Check if there's any site info in the map
                    if new_res is not None and new_pos is not None:
                        # Since we found the site in the site_map, add the
                        # updated site/position into the new site lists
                        # FIXME
                        if new_res == 'S':
                            new_mod_name = 'PhosphorylationSerine'
                        elif new_res == 'T':
                            new_mod_name = 'PhosphorylationThreonine'
                        elif new_res == 'Y':
                            new_mod_name = 'PhosphorylationTyrosine'
                        else:
                            raise Exception("Couldn't map residue.")
                        new_mod_list.append(new_mod_name)
                        new_modpos_list.append(new_pos)
                        # Add the mapped site to the invalid site list
                    invalid_sites[(agent.name, str(residue),
                                               str(old_modpos))] =  mapped_site
            # If the site is valid, add it and continue
            else:
                new_mod_list.append(old_mod)
                new_modpos_list.append(old_modpos)
        new_agent.mods = new_mod_list
        new_agent.mod_sites = new_modpos_list
        return (invalid_sites, new_agent)

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
