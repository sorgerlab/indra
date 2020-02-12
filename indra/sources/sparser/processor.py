from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import json
import logging
from copy import copy, deepcopy

from indra.util import read_unicode_csv
from indra.statements import *
from indra.databases import uniprot_client, hgnc_client

logger = logging.getLogger(__name__)


class SparserJSONProcessor(object):
    """Processor extracting INDRA Statements from Sparser's JSON output.

    Parameters
    ----------
    json_dict : list
        JSON containing the raw Sparser extractions.

    Attributes
    ----------
    json_stmts : list
        JSON containing the raw Sparser extractions.
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted by the processor.
    extraction_errors : list
        A list of tuples whose first element is the index of the JSON Statement
        in json_stmts and whose second element is an Exception that was
        thrown while processing the Statement with that index.
    """
    def __init__(self, json_dict):
        self.json_stmts = json_dict
        self.statements = []
        self.extraction_errors = []

    def get_statements(self):
        for idx, json_stmt in enumerate(self.json_stmts):
            try:
                # Step 1: Fix statement at the JSON level to avoid
                # deserialization issues
                json_stmt = fix_json_stmt(json_stmt)
                if json_stmt is None:
                    continue
                # Step 2: Deserialize into INDRA Statement
                stmt = Statement._from_json(json_stmt)
                # Step 3: Filter out invalid Statements
                if not check_statement_sanity(stmt):
                    continue
                # Step 4: Fix Agent names and grounding
                for agent in stmt.agent_list():
                    fix_agent(agent)
                # Step 5: Append to list of Statements
                self.statements.append(stmt)
            except SparserError as e:
                self.extraction_errors.append((idx, e))
            except NotAStatementName as e:
                self.extraction_errors.append((idx, e))
            except InvalidResidueError as e:
                self.extraction_errors.append((idx, e))
            except Exception as e:
                # Keep an eye on these and try to fix them as they come up, but
                # at least a reading job won't fail because of a couple
                # glitches. The logs should be processed, and processing errors
                # should be extracted and reported.
                logger.error("PROCESSING ERROR: Could not process json:\n%s"
                             % json.dumps(json_stmt))
                logger.exception(e)
                logger.error("END PROCESSING ERROR --------")
                self.extraction_errors.append((idx, e))

    def set_statements_pmid(self, pmid):
        """Set the evidence PMID of Statements that have been extracted.

        Parameters
        ----------
        pmid : str or None
            The PMID to be used in the Evidence objects of the Statements
            that were extracted by the processor.
        """
        # Replace PMID value in JSON dict first
        for stmt in self.json_stmts:
            evs = stmt.get('evidence', [])
            for ev in evs:
                ev['pmid'] = pmid
        # Replace PMID value in extracted Statements next
        for stmt in self.statements:
            for ev in stmt.evidence:
                ev.pmid = pmid

    def get_error_stats(self):
        """Return a Counter of error types that were encountered.

        Returns
        -------
        collections.Counter
            A Counter of error types that were encountered during processing.
        """
        from collections import Counter
        return Counter([e.__class__.__name__
                        for _, e in self.extraction_errors])


def fix_agent(agent):
    """Fix naming and grounding issues in an Agent, changes Agent in place."""
    if agent is None:
        return
    # First we fix some name spaces
    db_refs_tmp = copy(agent.db_refs)
    for db_ns, db_id in agent.db_refs.items():
        # Make sure CHEBI prefix is there
        if db_ns == 'CHEBI':
            if not db_id.startswith('CHEBI:'):
                db_refs_tmp['CHEBI'] = 'CHEBI:%s' % db_id
        # Change FA name space
        elif db_ns == 'FA':
            db_refs_tmp.pop('FA', None)
            db_refs_tmp['NXPFA'] = db_id
        # Change IPR name space
        elif db_ns == 'IPR':
            db_refs_tmp.pop('IPR', None)
            db_refs_tmp['IP'] = db_id
        # Change XFAM name space
        elif db_ns == 'XFAM':
            db_refs_tmp.pop('XFAM', None)
            db_refs_tmp['PF'] = db_id.split('.')[0]
        elif db_ns == 'GO':
            if db_id.startswith('GO:'):
                db_refs_tmp['GO'] = db_id
            else:
                db_refs_tmp['GO'] = 'GO:' + db_id
        # Change PCID name space
        elif db_ns == 'PCID':
            db_refs_tmp.pop('PCID', None)
            db_refs_tmp['PUBCHEM'] = db_id
    agent.db_refs = db_refs_tmp
    # Check if we have a FPLX entry and handle old BE mappings
    if 'BE' in agent.db_refs:
        agent.db_refs['FPLX'] = agent.db_refs.pop('BE')
    be_id = agent.db_refs.get('FPLX')
    # Try to map to FPLX from NXP, IPR, PF, NCIT
    if not be_id:
        for db_ns, db_id in agent.db_refs.items():
            be_id = famplex_map.get((db_ns, db_id))
            if be_id:
                break
    # Try mapping NCIT to specific genes if possible
    if not be_id and 'NCIT' in agent.db_refs:
        target = ncit_map.get(agent.db_refs['NCIT'])
        if target:
            agent.db_refs[target[0]] = target[1]
    # If the name is an UP ID, change it
    if agent.name and 'UP' not in agent.db_refs \
            and 'FPLX' not in agent.db_refs:
        if uniprot_client.get_gene_name(agent.name, web_fallback=False):
            agent.db_refs['UP'] = agent.name

    # Check what entries we have
    up_id = agent.db_refs.get('UP')
    hgnc_id = agent.db_refs.get('HGNC')
    # This is a special case that happens sometimes where agent.name is 'UP:
    # db_refs['UP'] is an empty string, and there is no other grounding.
    # In this case, we remove the empty UP grounding and reset the name to the
    # agent text.
    if not be_id and not hgnc_id and up_id == '':
        agent.name = agent.db_refs.get('TEXT', agent.name)
        agent.db_refs.pop('UP')
    # FPLX takes precedence if we have it
    elif be_id:
        agent.db_refs['FPLX'] = be_id
        agent.name = be_id
    elif hgnc_id:
        gene_name = hgnc_client.get_hgnc_name(hgnc_id)
        if gene_name:
            agent.name = gene_name
        if not up_id:
            up_id = hgnc_client.get_uniprot_id(hgnc_id)
            if up_id:
                if ', ' in up_id:
                    up_id = up_id.split(', ')[0]
                agent.db_refs['UP'] = up_id
    elif up_id:
        hgnc_id = uniprot_client.get_hgnc_id(up_id)
        if hgnc_id:
            agent.db_refs['HGNC'] = hgnc_id
            agent.name = hgnc_client.get_hgnc_name(hgnc_id)
        else:
            gene_name = uniprot_client.get_gene_name(up_id, web_fallback=False)
            if gene_name:
                agent.name = gene_name
            # If it doesn't have a gene name, it's better to just
            # use the raw string name otherwise Sparser sets
            # has Uniprot IDs or mnemonics as the name
            else:
                name = agent.db_refs.get('TEXT', agent.name)
                agent.name = name


def fix_json_stmt(json_stmt):
    # Step 1: fix JSON directly to reduce errors when deserializing
    # 1.1 Fix statement type issues
    # Change to IncreaseAmount
    if json_stmt['type'] == 'GeneTranscriptExpress':
        json_stmt['type'] = 'IncreaseAmount'
    stmt_type = json_stmt['type']
    stmt_class = get_statement_by_name(stmt_type)

    # 1.2 - Check for string agents.
    for ag_key in stmt_class._agent_order:
        json_stmt[ag_key] = fix_json_agent(json_stmt.get(ag_key))

    # 1.3 - Fix other misc things that are statement type specific
    if stmt_type in mod_class_names:
        position = json_stmt.get('position')
        residue = json_stmt.get('residue')
        if isinstance(position, list):
            if len(position) != 1:
                logger.error('Invalid position: %s' % position)
            else:
                json_stmt['position'] = position[0]
        if isinstance(residue, list):
            if len(residue) != 1:
                logger.error('Invalid residue: %s' % residue)
            elif not isinstance(residue[0], str):
                logger.error('Invalid residue: %s' % residue)
            else:
                json_stmt['residue'] = residue[0]
        elif isinstance(residue, bool):
            json_stmt['residue'] = None
    elif stmt_type in ('Activation', 'Inhibition'):
        obj_activity = json_stmt.get('obj_activity')
        if isinstance(obj_activity, list):
            if len(obj_activity) != 1:
                logger.error('Invalid object activity: %s'
                             % obj_activity)
            else:
                json_stmt['obj_activity'] = obj_activity[0]
        obj = json_stmt.get('obj')
        if isinstance(obj, (list, str)):
            raise InvalidAgent
    elif stmt_type == 'Translocation':
        # Fix locations if possible
        for loc_param in ('from_location', 'to_location'):
            loc = json_stmt.get(loc_param)
            if loc:
                # Some invalid locations are produced very often and can be
                # silently skipped
                if loc in known_invalid_locations:
                    loc = None
                # Otherwise we try to get a valid location via a mapping
                try:
                    loc = get_valid_location(loc)
                except InvalidLocationError:
                    logger.error('Invalid location: %s' % loc)
                    loc = None
                json_stmt[loc_param] = loc
        # Skip Translocation with both locations None
        if (json_stmt.get('from_location') is None
                and json_stmt.get('to_location') is None):
            raise TranslocationWithoutLocations
    elif stmt_type == 'IncreaseAmount':
        # Skip if there is no subject
        subj = json_stmt.get('subj')
        if not subj:
            raise MissingSubj

    # 1.4 - Fix evidence
    evs = json_stmt.get('evidence')
    if evs and isinstance(evs, list):
        ev = evs[0]
        text = ev.get('text')
        if not isinstance(text, str) or not text:
            raise MissingEvidenceText
        pmid = ev.get('pmid')
        if isinstance(pmid, str) and pmid.startswith('PMID'):
            json_stmt['evidence'][0]['pmid'] = \
                json_stmt['evidence'][0]['pmid'][4:]
    else:
        raise InvalidEvidence

    return json_stmt


def fix_json_agent(ag_obj):
    """Fix the JSON representation of an agent or list of agents."""
    if isinstance(ag_obj, str):
        logger.info("Fixing string agent: %s." % ag_obj)
        return {'name': ag_obj, 'db_refs': {'TEXT': ag_obj}}
    elif isinstance(ag_obj, list):
        # Recursive for complexes and similar.
        return [fix_json_agent(ag) for ag in ag_obj]
    elif isinstance(ag_obj, dict):
        # If an Agent is missing a name but has a TEXT db_refs then we set
        # that text as the name.
        if 'name' not in ag_obj:
            if 'TEXT' in ag_obj.get('db_refs', {}):
                ret = deepcopy(ag_obj)
                ret['name'] = ag_obj['db_refs']['TEXT']
                return ret
    return ag_obj


def check_statement_sanity(stmt):
    """Return True if the statement passes some sanity checks."""
    # Skip Statement if all agents are None
    if not any(stmt.agent_list()):
        raise NoAgents
    # Skip RegulateActivity if object is None
    if isinstance(stmt, RegulateActivity):
        if stmt.obj is None:
            raise MissingObj
        if stmt.subj is None:
            raise MissingSubj
    if isinstance(stmt, Modification):
        if stmt.sub is None:
            raise MissingSubj
    # Skip Complexes with less than 2 members
    if isinstance(stmt, Complex):
        if len(stmt.members) < 2:
            raise UnaryComplex
    return True


def _read_ncit_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../../resources/ncit_map.tsv')
    ncit_map = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    next(csv_rows)
    for row in csv_rows:
        ncit_id = row[0]
        target_ns = row[1]
        target_id = row[2]
        ncit_map[ncit_id] = (target_ns, target_id)
    return ncit_map


ncit_map = _read_ncit_map()


def _read_famplex_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../../resources/famplex_map.tsv')
    famplex_map = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        source_ns = row[0]
        source_id = row[1]
        be_id = row[2]
        if source_ns == 'NXP':
            source_ns = 'NXPFA'
            source_id = source_id.split(':')[1]
        famplex_map[(source_ns, source_id)] = be_id
    return famplex_map


famplex_map = _read_famplex_map()
mod_class_names = [cls.__name__ for cls in modclass_to_modtype.keys()]

known_invalid_locations = {'CHEBI:36080', 'GO:0038023'}


class SparserError(ValueError):
    pass


class InvalidAgent(SparserError):
    pass


class TranslocationWithoutLocations(SparserError):
    pass


class MissingSubj(SparserError):
    pass


class MissingObj(SparserError):
    pass


class UnaryComplex(SparserError):
    pass


class NoAgents(SparserError):
    pass


class MissingEvidenceText(SparserError):
    pass


class InvalidEvidence(SparserError):
    pass
