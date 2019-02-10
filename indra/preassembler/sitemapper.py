from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import os
import pickle
import logging
import textwrap
import requests
from copy import deepcopy
from functools import lru_cache
from protmapper import phosphosite_client
from protmapper.api import ProtMapper, default_site_map
from indra.statements import *
from indra.util import read_unicode_csv
from indra.config import has_config, get_config
from indra.databases import uniprot_client, hgnc_client
# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger(__name__)


class MappedStatement(object):
    """Information about a Statement found to have invalid sites.

    Parameters
    ----------
    original_stmt : :py:class:`indra.statements.Statement`
        The statement prior to mapping.
    mapped_mods : list of MappedSite
        A list of MappedSite objects.
    mapped_stmt : :py:class:`indra.statements.Statement`
        The statement after mapping. Note that if no information was found
        in the site map, it will be identical to the original statement.
    """
    def __init__(self, original_stmt, mapped_mods, mapped_stmt):
        self.original_stmt = original_stmt
        self.mapped_mods = mapped_mods
        self.mapped_stmt = mapped_stmt

    @python_2_unicode_compatible
    def __str__(self):
        if not self.mapped_mods:
            mm_str = str(self.mapped_mods)
        else:
            mm_ws = '\n' + (' ' * 17)
            mm_str = mm_ws.join([str(mm) for mm in self.mapped_mods])

        summary = textwrap.dedent("""
            MappedStatement:
                original_stmt: {0}
                mapped_mods: {1}
                mapped_stmt: {2}
            """)
        return summary.format(self.original_stmt, mm_str, self.mapped_stmt)

    def __repr__(self):
        return str(self)


class SiteMapper(ProtMapper):
    """
    Use curated site information to standardize modification sites in stmts.

    Parameters
    ----------
    site_map : dict (as returned by :py:func:`load_site_map`)
        A dict mapping tuples of the form `(gene, orig_res, orig_pos)` to a
        tuple of the form `(correct_res, correct_pos, comment)`, where `gene`
        is the string name of the gene (canonicalized to HGNC); `orig_res` and
        `orig_pos` are the residue and position to be mapped; `correct_res` and
        `correct_pos` are the corrected residue and position, and `comment` is
        a string describing the reason for the mapping (species error, isoform
        error, wrong residue name, etc.).
    use_cache : Optional[bool]
        If True, the SITEMAPPER_CACHE_PATH from the config (or environment)
        is loaded and cached mappings are read and written to the given path.
        Otherwise, no cache is used. Default: False

    Examples
    --------
    Fixing site errors on both the modification state of an agent (MAP2K1) and
    the target of a Phosphorylation statement (MAPK1):

    >>> map2k1_phos = Agent('MAP2K1', db_refs={'UP':'Q02750'}, mods=[
    ... ModCondition('phosphorylation', 'S', '217'),
    ... ModCondition('phosphorylation', 'S', '221')])
    >>> mapk1 = Agent('MAPK1', db_refs={'UP':'P28482'})
    >>> stmt = Phosphorylation(map2k1_phos, mapk1, 'T','183')
    >>> (valid, mapped) = default_mapper.map_sites([stmt])
    >>> valid
    []
    >>> mapped  # doctest:+IGNORE_UNICODE
    [
    MappedStatement:
        original_stmt: Phosphorylation(MAP2K1(mods: (phosphorylation, S, 217), (phosphorylation, S, 221)), MAPK1(), T, 183)
        mapped_mods: (('MAP2K1', 'S', '217'), ('S', '218', 'off by one'))
                     (('MAP2K1', 'S', '221'), ('S', '222', 'off by one'))
                     (('MAPK1', 'T', '183'), ('T', '185', 'off by two; mouse sequence'))
        mapped_stmt: Phosphorylation(MAP2K1(mods: (phosphorylation, S, 218), (phosphorylation, S, 222)), MAPK1(), T, 185)
    ]
    >>> ms = mapped[0]
    >>> ms.original_stmt
    Phosphorylation(MAP2K1(mods: (phosphorylation, S, 217), (phosphorylation, S, 221)), MAPK1(), T, 183)
    >>> ms.mapped_mods # doctest:+IGNORE_UNICODE
    [(('MAP2K1', 'S', '217'), ('S', '218', 'off by one')), (('MAP2K1', 'S', '221'), ('S', '222', 'off by one')), (('MAPK1', 'T', '183'), ('T', '185', 'off by two; mouse sequence'))]
    >>> ms.mapped_stmt
    Phosphorylation(MAP2K1(mods: (phosphorylation, S, 218), (phosphorylation, S, 222)), MAPK1(), T, 185)
    """
    def map_stmt_sites(self, stmt, do_methionine_offset=True,
                       do_orthology_mapping=True, do_isoform_mapping=True):
        stmt_copy = deepcopy(stmt)
        # For all statements, replace agents with invalid modifications
        mapped_sites = []
        new_agent_list = []

        for agent in stmt.agent_list():
            # If None agent, add to list and continue
            if agent is None:
                new_agent_list.append(agent)
                continue
            # Otherwise do mapping
            agent_mapped_sites, new_agent = self._map_agent_sites(
                agent,
                do_methionine_offset=do_methionine_offset,
                do_orthology_mapping=do_orthology_mapping,
                do_isoform_mapping=do_isoform_mapping
                )
            mapped_sites += agent_mapped_sites

            # Site map agents in the bound conditions
            for ind in range(len(new_agent.bound_conditions)):
                b = new_agent.bound_conditions[ind].agent
                agent_mapped_sites, new_b = self._map_agent_sites(
                    b,
                    do_methionine_offset=do_methionine_offset,
                    do_orthology_mapping=do_orthology_mapping,
                    do_isoform_mapping=do_isoform_mapping
                    )
                mapped_sites += agent_mapped_sites
                new_agent.bound_conditions[ind].agent = new_b

            new_agent_list.append(new_agent)

        # If we made any actual mappings, we set the new agent list which
        # includes any mapped agents
        if any([ms.has_mapping() for ms in mapped_sites]):
            stmt_copy.set_agent_list(new_agent_list)

        # --- Special handling for these statements ---
        # For modifications, fix residue and position
        if (isinstance(stmt, Modification) or
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
            old_mod = stmt._get_mod_condition()
            # Figure out if this site is invalid
            stmt_mapped_site = self._map_agent_mod(
                    agent_to_check,
                    old_mod,
                    do_methionine_offset=do_methionine_offset,
                    do_orthology_mapping=do_orthology_mapping,
                    do_isoform_mapping=do_isoform_mapping
                    )
            if stmt_mapped_site.has_mapping():
                stmt_copy.residue = stmt_mapped_site.mapped_res
                stmt_copy.position = stmt_mapped_site.mapped_pos
                mapped_sites.append(stmt_mapped_site)
        if any([ms is not None for ms in mapped_sites]):
            mapped_stmt = MappedStatement(stmt, mapped_sites, stmt_copy)
        else:
            mapped_stmt = None
        return mapped_stmt

    def map_sites(self, stmts, do_methionine_offset=True,
                  do_orthology_mapping=True, do_isoform_mapping=True):
        """Check a set of statements for invalid modification sites.

        Statements are checked against Uniprot reference sequences to determine
        if residues referred to by post-translational modifications exist at
        the given positions.

        If there is nothing amiss with a statement (modifications on any of the
        agents, modifications made in the statement, etc.), then the statement
        goes into the list of valid statements. If there is a problem with the
        statement, the offending modifications are looked up in the site map
        (:py:attr:`site_map`), and an instance of :py:class:`MappedStatement`
        is added to the list of mapped statements.

        Parameters
        ----------
        stmts : list of :py:class:`indra.statement.Statement`
            The statements to check for site errors.
        do_methionine_offset : boolean
            Whether to check for off-by-one errors in site position (possibly)
            attributable to site numbering from mature proteins after
            cleavage of the initial methionine. If True, checks the reference
            sequence for a known modification at 1 site position greater
            than the given one; if there exists such a site, creates the
            mapping. Default is True.
        do_orthology_mapping : boolean
            Whether to check sequence positions for known modification sites
            in mouse or rat sequences (based on PhosphoSitePlus data). If a
            mouse/rat site is found that is linked to a site in the human
            reference sequence, a mapping is created. Default is True.
        do_isoform_mapping : boolean
            Whether to check sequence positions for known modifications
            in other human isoforms of the protein (based on PhosphoSitePlus
            data). If a site is found that is linked to a site in the human
            reference sequence, a mapping is created. Default is True.

        Returns
        -------
        tuple
            2-tuple containing (valid_statements, mapped_statements). The first
            element of the tuple is a list valid statements
            (:py:class:`indra.statement.Statement`) that were not found to
            contain any site errors. The second element of the tuple is a list
            of mapped statements (:py:class:`MappedStatement`) with information
            on the incorrect sites and corresponding statements with correctly
            mapped sites.
        """
        valid_statements = []
        mapped_statements = []

        for stmt in stmts:
            mapped_stmt = self.map_stmt_sites(stmt, do_methionine_offset,
                                              do_orthology_mapping,
                                              do_isoform_mapping)

            # If the invalid_sites list isn't empty, that means that there were
            # incorrect residues for this statement; add the statement to
            # the mapped_statements list
            if mapped_stmt is not None:
                mapped_statements.append(mapped_stmt)
            else:
                valid_statements.append(stmt)

        return valid_statements, mapped_statements

    def _map_agent_sites(self, agent, do_methionine_offset=True,
                         do_orthology_mapping=True, do_isoform_mapping=True):
        """Check an agent for invalid sites and update if necessary.

        Parameters
        ----------
        agent : :py:class:`indra.statements.Agent`
            Agent to check for invalid modification sites.
        do_methionine_offset : boolean
            Whether to check for off-by-one errors in site position (possibly)
            attributable to site numbering from mature proteins after
            cleavage of the initial methionine. If True, checks the reference
            sequence for a known modification at 1 site position greater
            than the given one; if there exists such a site, creates the
            mapping. Default is True.
        do_orthology_mapping : boolean
            Whether to check sequence positions for known modification sites
            in mouse or rat sequences (based on PhosphoSitePlus data). If a
            mouse/rat site is found that is linked to a site in the human
            reference sequence, a mapping is created. Default is True.
        do_isoform_mapping : boolean
            Whether to check sequence positions for known modifications
            in other human isoforms of the protein (based on PhosphoSitePlus
            data). If a site is found that is linked to a site in the human
            reference sequence, a mapping is created. Default is True.

        Returns
        -------
        tuple
            The first element is a list of MappedSite objects, the second
            element is either the original Agent, if unchanged, or a copy
            of it.
        """
        if agent is None:
            return [], agent
        new_agent = deepcopy(agent)
        # If there are no modifications on this agent, then we can return the
        # copy of the agent
        if not agent.mods:
            return [], new_agent
        mapped_sites = []
        for idx, mod_condition in enumerate(agent.mods):
            mapped_site = \
                self._map_agent_mod(agent, mod_condition,
                                    do_methionine_offset=do_methionine_offset,
                                    do_orthology_mapping=do_orthology_mapping,
                                    do_isoform_mapping=do_isoform_mapping)
            # If we couldn't do the mapping or the mapped site isn't invalid
            # then we don't need to change the existing ModCondition
            if not mapped_site or mapped_site.not_invalid():
                continue
            # Otherwise, if there is a mapping, we replace the old ModCondition
            # with the new one where only the residue and position are updated,
            # the mod type and the is modified flag are kept.
            if mapped_site.has_mapping():
                mc = ModCondition(mod_condition.mod_type,
                                  mapped_site.mapped_res,
                                  mapped_site.mapped_pos,
                                  mod_condition.is_modified)
                new_agent.mods[idx] = mc
            # Finally, whether or not we have a mapping, we keep track of mapped
            # sites and make them available to the caller
            mapped_sites.append(mapped_site)
        return mapped_sites, new_agent

    def _map_agent_mod(self, agent, mod_condition, do_methionine_offset=True,
                       do_orthology_mapping=True, do_isoform_mapping=True):
        """Check a single modification site on an agent and look for a mapping.

        Parameters
        ----------
        agent : :py:class:`indra.statements.Agent`
            Agent to check for invalid modification sites.
        mod_condition : :py:class:`indra.statements.ModCondition`
            Modification to check for validity and map.
        do_methionine_offset : boolean
            Whether to check for off-by-one errors in site position (possibly)
            attributable to site numbering from mature proteins after
            cleavage of the initial methionine. If True, checks the reference
            sequence for a known modification at 1 site position greater
            than the given one; if there exists such a site, creates the
            mapping. Default is True.
        do_orthology_mapping : boolean
            Whether to check sequence positions for known modification sites
            in mouse or rat sequences (based on PhosphoSitePlus data). If a
            mouse/rat site is found that is linked to a site in the human
            reference sequence, a mapping is created. Default is True.
        do_isoform_mapping : boolean
            Whether to check sequence positions for known modifications
            in other human isoforms of the protein (based on PhosphoSitePlus
            data). If a site is found that is linked to a site in the human
            reference sequence, a mapping is created. Default is True.

        Returns
        -------
        protmapper.MappedSite or None
            A MappedSite object.
        """
        # Get the UniProt ID of the agent, if not found, return
        up_id = _get_uniprot_id(agent)
        if not up_id:
            logger.debug("No uniprot ID for %s" % agent.name)
            return None
        # If no site information for this residue, skip
        if mod_condition.position is None or mod_condition.residue is None:
            return None
        # Otherwise, try to map it and return the mapped site
        mapped_site = \
            self.map_to_human_ref(up_id, 'uniprot',
                                  mod_condition.residue,
                                  mod_condition.position,
                                  do_methionine_offset=do_methionine_offset,
                                  do_orthology_mapping=do_orthology_mapping,
                                  do_isoform_mapping=do_isoform_mapping)
        return mapped_site


default_mapper = SiteMapper(default_site_map)


@lru_cache(maxsize=10000)
def _get_uniprot_id(agent):
    """Get the Uniprot ID for an agent, looking up in HGNC if necessary.

    If the Uniprot ID is a list then return the first ID by default.
    """
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
