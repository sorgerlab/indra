__all__ = ['GroundingMapper', 'load_grounding_map', 'default_grounding_map',
           'default_agent_map', 'default_ignores', 'default_misgrounding_map',
           'default_mapper', 'gm']
import os
import csv
import json
import logging
from copy import deepcopy
from indra.statements import Agent
from indra.databases import hgnc_client
from indra.util import read_unicode_csv
from indra.preassembler.grounding_mapper.gilda import get_gilda_models
from indra.ontology.standardize import standardize_db_refs, \
    standardize_agent_name
from .disambiguate import adeft_disambiguators, DisambManager

logger = logging.getLogger(__name__)


class GroundingMapper(object):
    """Maps grounding of INDRA Agents based on a given grounding map.

    Each parameter, if not provided will result in loading the corresponding
    built-in grounding resource. To explicitly avoid loading the default,
    pass in an empty data structure as the given parameter, e.g., ignores=[].

    Parameters
    ----------
    grounding_map : Optional[dict]
        The grounding map, a dictionary mapping strings (entity names) to
        a dictionary of database identifiers.
    agent_map : Optional[dict]
        A dictionary mapping strings to grounded INDRA Agents with given state.
    ignores : Optional[list]
        A list of entity strings that, if encountered will result in the
        corresponding Statement being discarded.
    misgrounding_map : Optional[dict]
        A mapping dict similar to the grounding map which maps entity strings
        to a given grounding which is known to be incorrect and should be
        removed if encountered (making the remaining Agent ungrounded).
    use_adeft : Optional[bool]
        If True, Adeft will be attempted to be used for disambiguation of
        acronyms. Default: True
    gilda_mode : Optional[str]
        If None, Gilda will not be used at all. If 'web', the GILDA_URL
        setting from the config file or as an environmental variable
        is assumed to be the web service endpoint through which Gilda is used.
        If 'local', we assume that the gilda Python package is installed
        and will be used.
    """
    def __init__(self, grounding_map=None, agent_map=None, ignores=None,
                 misgrounding_map=None, use_adeft=True, gilda_mode=None):
        self.grounding_map = grounding_map if grounding_map is not None \
            else default_grounding_map
        self.check_grounding_map(self.grounding_map)
        self.agent_map = agent_map if agent_map is not None \
            else default_agent_map
        self.ignores = set(ignores) if ignores else default_ignores
        self.misgrounding_map = misgrounding_map if misgrounding_map \
            else default_misgrounding_map
        self.use_adeft = use_adeft
        self.disamb_manager = DisambManager()
        self.gilda_mode = gilda_mode
        self._gilda_models = None

    @property
    def gilda_models(self):
        if self._gilda_models is None:
            self._gilda_models = get_gilda_models(self.gilda_mode) \
                if self.gilda_mode else []
        return self._gilda_models

    @gilda_models.setter
    def gilda_models(self, models):
        self._gilda_models = models

    @staticmethod
    def check_grounding_map(gm):
        """Run sanity checks on the grounding map, raise error if needed."""
        for key, refs in gm.items():
            if not refs:
                continue
            if 'HGNC' in refs and \
                    hgnc_client.get_hgnc_name(refs['HGNC']) is None:
                raise ValueError('HGNC:%s for key %s in the grounding map is '
                                 'not a valid ID' % (refs['HGNC'], key))

    def map_stmts(self, stmts, do_rename=True):
        """Return a new list of statements whose agents have been mapped

        Parameters
        ----------
        stmts : list of :py:class:`indra.statements.Statement`
            The statements whose agents need mapping
        do_rename: Optional[bool]
            If True, the Agent name is updated based on the mapped grounding.
            If do_rename is True the priority for setting the name is
            FamPlex ID, HGNC symbol, then the gene name
            from Uniprot. Default: True

        Returns
        -------
        mapped_stmts : list of :py:class:`indra.statements.Statement`
            A list of statements given by mapping the agents from each
            statement in the input list
        """
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

    def map_agents_for_stmt(self, stmt, do_rename=True):
        """Return a new Statement whose agents have been grounding mapped.

        Parameters
        ----------
        stmt : :py:class:`indra.statements.Statement`
            The Statement whose agents need mapping.
        do_rename: Optional[bool]
            If True, the Agent name is updated based on the mapped grounding.
            If do_rename is True the priority for setting the name is
            FamPlex ID, HGNC symbol, then the gene name
            from Uniprot. Default: True

        Returns
        -------
        mapped_stmt : :py:class:`indra.statements.Statement`
            The mapped Statement.
        """
        mapped_stmt = deepcopy(stmt)

        # Iterate over the agents
        # Update agents directly participating in the statement
        agent_list = mapped_stmt.agent_list()
        for idx, agent in enumerate(agent_list):
            # If the agent is None, we do nothing
            if agent is None:
                continue
            # If the agent's TEXT is in the ignores list, we return None to
            # then filter out the Statement
            agent_txts = {agent.db_refs[t] for t in {'TEXT', 'TEXT_NORM'}
                          if t in agent.db_refs}
            if agent_txts and agent_txts & set(self.ignores):
                return None

            # Check if an adeft model exists for agent text
            adeft_success = False
            if self.use_adeft and agent_txts and agent_txts & \
                    set(adeft_disambiguators):
                try:
                    # Us the longest match for disambiguation
                    txt_for_adeft = sorted(agent_txts &
                                           set(adeft_disambiguators),
                                           key=lambda x: len(x))[-1]
                    adeft_success = self.disamb_manager.\
                        run_adeft_disambiguation(mapped_stmt, agent, idx,
                                                 txt_for_adeft)
                except Exception as e:
                    logger.error('There was an error during Adeft'
                                 ' disambiguation of %s.' % str(agent_txts))
                    logger.error(e)

            gilda_success = False
            # Gilda is not used if agent text is in the grounding map
            if not adeft_success and self.gilda_mode and \
               not agent_txts & set(self.grounding_map) and \
               agent_txts & set(self.gilda_models):
                try:
                    # Us the longest match for disambiguation
                    txt_for_gilda = sorted(agent_txts & set(self.gilda_models),
                                           key=lambda x: len(x))[-1]
                    gilda_success = self.disamb_manager.\
                        run_gilda_disambiguation(mapped_stmt, agent, idx,
                                                 txt_for_gilda,
                                                 mode=self.gilda_mode)
                except Exception as e:
                    logger.error('There was an error during Gilda'
                                 ' disambiguation of %s.' % str(agent_txts))
                    logger.error(e)

            # If Adeft and Gilda were not used or didn't succeed, we do
            # grounding mapping
            new_agent = self.map_agent(agent, do_rename) \
                if not (adeft_success or gilda_success) else agent

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
                    bc.agent = self.map_agent(bc.agent, do_rename)
                    if not bc.agent:
                        # Skip the entire statement if the agent maps to None
                        # in the grounding map
                        return None

        return mapped_stmt

    def map_agent(self, agent, do_rename):
        """Return the given Agent with its grounding mapped.

        This function grounds a single agent. It returns the new Agent object
        (which might be a different object if we load a new agent state
        from json) or the same object otherwise.

        Parameters
        ----------
        agent : :py:class:`indra.statements.Agent`
            The Agent to map.
        do_rename: bool
            If True, the Agent name is updated based on the mapped grounding.
            If do_rename is True the priority for setting the name is
            FamPlex ID, HGNC symbol, then the gene name
            from Uniprot.

        Returns
        -------
        grounded_agent : :py:class:`indra.statements.Agent`
            The grounded Agent.
        """
        # We always standardize DB refs as a functionality in the
        # GroundingMapper. If a new module is implemented which is
        # responsible for standardizing grounding, this can be removed.
        agent.db_refs = self.standardize_db_refs(agent.db_refs)
        # If there is no TEXT available, we can return immediately since we
        # can't do mapping
        agent_txts = sorted({agent.db_refs[t] for t in {'TEXT', 'TEXT_NORM'}
                             # Note that get here will correctly handle both
                             # a non-existent entry and a None entry which
                             # sometimes appears
                             if agent.db_refs.get(t)}, key=lambda x: len(x),
                            reverse=True)
        if not agent_txts:
            # We still do the name standardization here
            if do_rename:
                self.standardize_agent_name(agent, standardize_refs=False)
            return agent

        # 1. Check if there is a full agent mapping and apply if there is
        for agent_text in agent_txts:
            if agent_text in self.agent_map:
                mapped_to_agent = \
                    Agent._from_json(self.agent_map[agent_text]['agent'])
                return mapped_to_agent

        # 2. Look agent text up in the grounding map
        for agent_text in agent_txts:
            if agent_text in self.grounding_map:
                self.update_agent_db_refs(agent, self.grounding_map[agent_text],
                                          do_rename)
                return agent

        # 3. Look agent text up in the misgrounding map
        for agent_text in agent_txts:
            if agent_text in self.misgrounding_map:
                self.remove_agent_db_refs(agent,
                                          self.misgrounding_map[agent_text])
        # This happens when there is an Agent text but it is not in the
        # grounding map. We still do the name standardization here.
        if do_rename:
            self.standardize_agent_name(agent, standardize_refs=False)
        # Otherwise just return
        return agent

    def update_agent_db_refs(self, agent, db_refs, do_rename=True):
        """Update db_refs of agent using the grounding map

        If the grounding map is missing one of the HGNC symbol or Uniprot ID,
        attempts to reconstruct one from the other.

        Parameters
        ----------
        agent : :py:class:`indra.statements.Agent`
            The agent whose db_refs will be updated
        db_refs : dict
            The db_refs so set for the agent.
        do_rename: Optional[bool]
            If True, the Agent name is updated based on the mapped grounding.
            If do_rename is True the priority for setting the name is
            FamPlex ID, HGNC symbol, then the gene name
            from Uniprot. Default: True
        """
        # Standardize the IDs in the db_refs dict and set it as the Agent's
        # db_refs
        txt = agent.db_refs.get('TEXT')
        agent.db_refs = self.standardize_db_refs(deepcopy(db_refs))
        if txt:
            agent.db_refs['TEXT'] = txt
        # Finally, if renaming is needed we standardize the Agent's name
        if do_rename:
            self.standardize_agent_name(agent, standardize_refs=False)

    def remove_agent_db_refs(self, agent, db_refs):
        # Standardize the IDs in the db_refs dict and set it as the Agent's
        # db_refs
        standard_refs = self.standardize_db_refs(deepcopy(db_refs))
        # If there is any overlap between the Agent's db_refs and the db_refs
        # that are to be eliminated, we consider the Agent's db_refs to be
        # invalid and remove them. We then reset the Agent's name to
        # its TEXT value if available.
        preserve_refs = {k: agent.db_refs[k] for k in {'TEXT', 'TEXT_NORM'}
                         if k in agent.db_refs}
        if set(standard_refs.items()) & set(agent.db_refs.items()):
            agent.db_refs = preserve_refs
            if 'TEXT_NORM' in agent.db_refs:
                agent.name = agent.db_refs['TEXT_NORM']
            elif 'TEXT' in agent.db_refs:
                agent.name = agent.db_refs['TEXT']

    @staticmethod
    def standardize_db_refs(db_refs):
        """Return a standardized db refs dict for a given db refs dict.

        Parameters
        ----------
        db_refs : dict
            A dict of db refs that may not be standardized, i.e., may be
            missing an available UP ID corresponding to an existing HGNC ID.

        Returns
        -------
        dict
            The db_refs dict with standardized entries.
        """
        return standardize_db_refs(db_refs)

    @staticmethod
    def standardize_agent_name(agent, standardize_refs=True):
        """Standardize the name of an Agent based on grounding information.

        If an agent contains a FamPlex grounding, the FamPlex ID is used as a
        name. Otherwise if it contains a Uniprot ID, an attempt is made to find
        the associated HGNC gene name. If one can be found it is used as the
        agent name and the associated HGNC ID is added as an entry to the
        db_refs. Similarly, CHEBI, MESH and GO IDs are used in this order of
        priority to assign a standardized name to the Agent. If no relevant
        IDs are found, the name is not changed.

        Parameters
        ----------
        agent : indra.statements.Agent
            An INDRA Agent whose name attribute should be standardized based
            on grounding information.
        standardize_refs : Optional[bool]
            If True, this function assumes that the Agent's db_refs need to
            be standardized, e.g., HGNC mapped to UP.
            Default: True
        """
        return standardize_agent_name(agent,
                                      standardize_refs=standardize_refs)

    @staticmethod
    def rename_agents(stmts):
        """Return a list of mapped statements with updated agent names.

        Creates a new list of statements without modifying the original list.

        Parameters
        ----------
        stmts : list of :py:class:`indra.statements.Statement`
            List of statements whose Agents need their names updated.

        Returns
        -------
        mapped_stmts : list of :py:class:`indra.statements.Statement`
            A new list of Statements with updated Agent names
        """
        # Make a copy of the stmts
        mapped_stmts = deepcopy(stmts)
        # Iterate over the statements
        for _, stmt in enumerate(mapped_stmts):
            # Iterate over the agents
            for agent in stmt.agent_list():
                GroundingMapper.standardize_agent_name(agent, True)
        return mapped_stmts


# TODO: handle the cases when there is more than one entry for the same
# key (e.g., ROS, ER)
def load_grounding_map(grounding_map_path, lineterminator='\r\n',
                       hgnc_symbols=True):
    """Return a grounding map dictionary loaded from a csv file.

    In the file pointed to by grounding_map_path, the number of name_space ID
    pairs can vary per row and commas are
    used to pad out entries containing fewer than the maximum amount of
    name spaces appearing in the file. Lines should be terminated with \r\n
    both a carriage return and a new line by default.

    Optionally, one can specify another csv file (pointed to by ignore_path)
    containing agent texts that are degenerate and should be filtered out.

    It is important to note that this function assumes that the mapping file
    entries for the HGNC key are symbols not IDs. These symbols are converted
    to IDs upon loading here.

    Parameters
    ----------
    grounding_map_path : str
        Path to csv file containing grounding map information. Rows of the file
        should be of the form <agent_text>,<name_space_1>,<ID_1>,...
        <name_space_n>,<ID_n>
    lineterminator : Optional[str]
        Line terminator used in input csv file. Default: \r\n
    hgnc_symbols : Optional[bool]
        Set to True if the grounding map file contains HGNC symbols rather than
        IDs. In this case, the entries are replaced by IDs. Default: True

    Returns
    -------
    g_map : dict
        The grounding map constructed from the given files.
    """
    gmap = {}
    map_rows = read_unicode_csv(grounding_map_path, delimiter=',',
                                quotechar='"',
                                quoting=csv.QUOTE_MINIMAL,
                                lineterminator=lineterminator)
    for row in map_rows:
        txt = row[0]
        keys = [entry for entry in row[1::2] if entry]
        values = [entry for entry in row[2::2] if entry]
        if not keys or not values:
            logger.warning('Missing grounding entries for %s, skipping.' % txt)
            continue
        if len(keys) != len(values):
            logger.warning('Mismatched keys and values in row %s, skipping.' %
                           str(row))
            continue
        gmap[txt] = dict(zip(keys, values))
    if hgnc_symbols:
        gmap = replace_hgnc_symbols(gmap)
    return gmap


def replace_hgnc_symbols(gmap):
    """Replace HGNC symbols with IDs in a grounding map."""
    for txt, mapped_refs in deepcopy(gmap).items():
        hgnc_sym = mapped_refs.get('HGNC')
        if hgnc_sym:
            hgnc_id = hgnc_client.get_hgnc_id(hgnc_sym)
            # Override the HGNC symbol entry from the grounding
            # map with an HGNC ID
            if hgnc_id:
                mapped_refs['HGNC'] = hgnc_id
            else:
                logger.error('No HGNC ID corresponding to gene '
                             'symbol %s in grounding map.' % hgnc_sym)
                # Remove the HGNC symbol in this case
                mapped_refs.pop('HGNC')
        # In case the only grounding was eliminated, we remove the entry
        # completely
        if mapped_refs:
            gmap[txt] = mapped_refs
    return gmap


def _get_resource_path(*suffixes):
    return os.path.join(os.path.dirname(__file__), os.pardir, os.pardir,
                        'resources', *suffixes)


def _load_default_grounding_map():
    default_grounding_map_path = \
        _get_resource_path('grounding', 'grounding_map.csv')
    gmap = load_grounding_map(default_grounding_map_path, hgnc_symbols=True)
    return gmap


def _load_default_misgrounding_map():
    default_misgrounding_map_path = \
        _get_resource_path('grounding', 'misgrounding_map.csv')
    gmap = load_grounding_map(default_misgrounding_map_path, hgnc_symbols=False)
    return gmap


def _load_default_agent_map():
    default_agent_grounding_path = \
        _get_resource_path('grounding', 'agents.json')
    with open(default_agent_grounding_path, 'r') as fh:
        agent_map = json.load(fh)
    return agent_map


def _load_default_ignores():
    default_ignore_path = _get_resource_path('grounding', 'ignore.csv')
    with open(default_ignore_path, 'r') as fh:
        ignores = [l.strip() for l in fh.readlines()]
    return ignores


default_grounding_map = _load_default_grounding_map()
gm = default_grounding_map  # For backwards compatibility, redundant
default_misgrounding_map = _load_default_misgrounding_map()
default_agent_map = _load_default_agent_map()
default_ignores = _load_default_ignores()
default_mapper = GroundingMapper(default_grounding_map,
                                 agent_map=default_agent_map,
                                 ignores=default_ignores,
                                 misgrounding_map=default_misgrounding_map)
