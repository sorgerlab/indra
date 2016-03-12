import os
import csv
from collections import namedtuple
from indra.pysb_assembler import abbrevs as site_abbrevs
from indra.databases import uniprot_client

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

        if isinstance(stmt, Complex):
            for m in stmt.members:
                failures += check_agent_mod(m)
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

    def check_agent_mod(self, agent, mods=None, mod_sites=None):
        """Look up the modification site in Uniprot and then the site map.
        """

        invalid_sites = {}
        # If no UniProt ID is found, we don't report a failure
        up_id = agent.db_refs.get('UP')
        if up_id is None:
            return invalid_sites

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
