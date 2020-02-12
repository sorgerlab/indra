from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import logging
import collections
from copy import copy, deepcopy

from indra.util import read_unicode_csv
from indra.literature import id_lookup
from indra.statements import *
from indra.databases import uniprot_client, hgnc_client

logger = logging.getLogger(__name__)


try:
    basestring
except:
    basestring = str


def _fix_json_agents(ag_obj):
    """Fix the json representation of an agent."""
    if isinstance(ag_obj, str):
        logger.info("Fixing string agent: %s." % ag_obj)
        ret = {'name': ag_obj, 'db_refs': {'TEXT': ag_obj}}
    elif isinstance(ag_obj, list):
        # Recursive for complexes and similar.
        ret = [_fix_json_agents(ag) for ag in ag_obj]
    elif isinstance(ag_obj, dict) and 'TEXT' in ag_obj.keys():
        ret = deepcopy(ag_obj)
        text = ret.pop('TEXT')
        ret['db_refs']['TEXT'] = text
    else:
        ret = ag_obj
    return ret


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
    """
    def __init__(self, json_dict):
        self.json_stmts = json_dict
        self.statements = []

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
            except NotAStatementName:
                logger.error("%s is not a valid Statement type." %
                             json_stmt.get('type'))
            '''
            except Exception as e:
                # Keep an eye on these and try to fix them as they come up, but
                # at least a reading job won't fail because of a couple
                # glitches. The logs should be processed, and processing errors
                # should be extracted and reported.
                logger.error("PROCESSING ERROR: Could not process json:\n%s"
                             % str(json_stmt))
                logger.exception(e)
                logger.error("END PROCESSING ERROR --------")
            '''

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


def fix_agent(agent):
    """Fix naming and grounding issues in an Agent, changes Agent in place."""
    if agent is None:
        return
    # First we fix some name spaces
    db_refs_tmp = copy(agent.db_refs)
    for db_ns, db_id in agent.db_refs.items():
        # Change FA name space
        if db_ns == 'FA':
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
        if uniprot_client.get_gene_name(agent.name):
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
                agent.db_refs['UP'] = up_id
    elif up_id:
        hgnc_id = uniprot_client.get_hgnc_id(up_id)
        if hgnc_id:
            agent.db_refs['HGNC'] = hgnc_id
            agent.name = hgnc_client.get_hgnc_name(hgnc_id)
        else:
            gene_name = uniprot_client.get_gene_name(up_id)
            if gene_name:
                agent.name = gene_name
            # If it doesn't have a gene name, it's better to just
            # use the raw string name otherwise Sparser sets
            # has Uniprot IDs or mnemonics as the name
            else:
                name = agent.db_refs.get('TEXT', agent.name)
                agent.name = name


def fix_json_stmt(json_stmt):
    # If type isn't there, we have a very serious problem.
    stmt_type = json_stmt['type']

    # Step 1: fix JSON directly to reduce errors when deserializing
    # 1.1 - Check for string agents.
    stmt_class = get_statement_by_name(stmt_type)
    for ag_key in stmt_class._agent_order:
        json_stmt[ag_key] = _fix_json_agents(json_stmt.get(ag_key))

    # 1.2 - Fix other misc things.
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
            elif not isinstance(residue[0], basestring):
                logger.error('Invalid residue: %s' % residue)
            else:
                json_stmt['residue'] = residue[0]
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
            return None
    elif stmt_type == 'Translocation':
        # Fix locations if possible
        for loc_param in ('from_location', 'to_location'):
            loc = json_stmt.get(loc_param)
            if loc:
                try:
                    loc = get_valid_location(loc)
                except InvalidLocationError:
                    logger.error('Invalid location: %s' % loc)
                    loc = None
                json_stmt[loc_param] = loc
        # Skip Translocation with both locations None
        if (json_stmt.get('from_location') is None
                and json_stmt.get('to_location') is None):
            return None
    elif stmt_type == 'GeneTranscriptExpress':
        # Skip if there is no subject
        subj = json_stmt.get('subj')
        if not subj:
            return None
        # Change to IncreaseAmount
        json_stmt['type'] = 'IncreaseAmount'
    return json_stmt


def check_statement_sanity(stmt):
    """Return True if the statement passes some sanity checks."""
    # Skip Statement if all agents are None
    if not any(stmt.agent_list()):
        return False
    # Skip RegulateActivity if object is None
    if isinstance(stmt, RegulateActivity):
        if stmt.obj is None or stmt.subj is None:
            return False
    if isinstance(stmt, Modification):
        if stmt.sub is None:
            return False
    # Skip Complexes with less than 2 members
    if isinstance(stmt, Complex):
        if len(stmt.members) < 2:
            return False
    return True


class SparserXMLProcessor(object):
    """Processor extracting INDRA Statements from Sparser's XML output.

    Parameters
    ----------
    xml_etree : xml.etree.ElementTree
        An ElementTree containing the XML output of Sparser.

    Attributes
    ----------
    tree : objectpath.Tree
        The objectpath Tree object representing the extractions.
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted by the processor.
    pmid : str
        The pmid of the content the statements were extracted from.
    """
    def __init__(self, xml_etree):
        self.tree = xml_etree
        self.statements = []
        # Extract all sems by category
        self._sems = collections.defaultdict(list)
        for interp in self.tree.findall('interpretation'):
            sentence = interp.find('sentence-text').text
            sems = interp.findall('sem')
            for sem in sems:
                ref = sem.find('ref')
                if ref is not None:
                    category = ref.attrib['category']
                    self._sems[category].append((sem, sentence))
        # Get citation info
        pmcid = self.tree.attrib.get('pmcid')
        pmid = self.tree.attrib.get('pmid')
        if not pmid:
            pmid = self.tree.attrib.get('id')
        self.pmid = None
        if pmid:
            if pmid.startswith('PMID'):
                pmid = pmid[4:]
            self.pmid = pmid
        elif pmcid:
            ids = id_lookup(pmcid, 'pmcid')
            pmid = ids.get('pmid')
            if pmid is not None:
                self.pmid = pmid

    def get_modifications(self):
        mod_events = {}
        mod_types = list(modtype_to_modclass.keys()) + \
            ['phosphorylate', 'dephosphorylate']
        for mod_type in mod_types:
            mod_events[mod_type] = self._sems.get(mod_type)
        if not mod_events:
            return
        for mod_type, sems in mod_events.items():
            if not sems:
                continue
            if mod_type in ['phosphorylate', 'dephosphorylate']:
                indra_class = modtype_to_modclass.get(mod_type[:-1] + 'ion')
            else:
                indra_class = modtype_to_modclass.get(mod_type)
            if indra_class is None:
                logger.warning('Unhandled modification type %s' % mod_type)
                continue
            for event, sentence in sems:
                # Get enzyme agent
                enzyme = event.find("ref/var/[@name='agent']/ref")
                if enzyme is None:
                    enz = None
                else:
                    enz = self._get_agent_from_ref(enzyme)

                # Get substrate agent
                substrate = event.find("ref/var/[@name='substrate']/ref")
                if substrate is None:
                    substrate = event.find("ref/var/" +
                                           "[@name='agent-or-substrate']/ref")
                    if substrate is None:
                        logger.debug('Skipping phosphorylation without' +
                                     'substrate.')
                        continue
                sub = self._get_agent_from_ref(substrate)
                if sub is None:
                    logger.debug('Skipping phosphorylation without substrate.')
                    continue

                # Get site
                residue = None
                position = None
                site = event.find("ref/var/[@name='site']")
                if site is not None:
                    residue, position = self._get_site(site)

                # Get evidence
                ev = self._get_evidence(sentence)

                st = indra_class(enz, sub, residue, position, evidence=[ev])
                self.statements.append(st)

    def get_activations(self):
        act_events = self._sems.get('bio-activate')
        inh_events = self._sems.get('inhibit')
        for event_type, events in zip(['act', 'inh'],
                                      [act_events, inh_events]):
            if not events:
                continue
            for event, sentence in events:
                # Get the subject of the activation
                subj_ref = event.find("ref/var/[@name='by-means-of-or-agent']"
                                      "/ref")
                if subj_ref is None:
                    subj_ref = event.find("ref/var/[@name='agent']/ref")
                if subj_ref is None:
                    logger.debug('Skipping activation without subject.')
                    continue
                subj = self._get_agent_from_ref(subj_ref)
                if subj is None:
                    logger.debug('Skipping activation without subject.')
                    continue
                # Get the object of the activation
                obj_ref = event.find("ref/var/[@name='object']/ref")
                if obj_ref is None:
                    logger.debug('Skipping activation without object.')
                    continue
                obj = self._get_agent_from_ref(obj_ref)
                if obj is None:
                    logger.debug('Skipping activation without object.')
                    continue
                # Get evidence
                ev = self._get_evidence(sentence)
                if event_type == 'act':
                    st = Activation(subj, obj, evidence=[ev])
                else:
                    st = Inhibition(subj, obj, evidence=[ev])
                self.statements.append(st)

    def get_translocations(self):
        pass

    def _get_agent_from_ref(self, ref):
        # TODO: handle collections
        if ref.attrib.get('category') == 'collection':
            #logger.warning('Skipping collection Agent.')
            return None

        # Find the name, uid and raw-text tags first and get their text
        # content if available
        uid_tag = ref.find("var/[@name='uid']")
        name_tag = ref.find("var/[@name='name']")
        text_tag = ref.find("var/[@name='raw-text']")
        if name_tag is not None and name_tag.text:
            name = name_tag.text
        else:
            name = None
        if uid_tag is not None and uid_tag.text:
            uid = uid_tag.text
        else:
            uid = None
        if text_tag is not None and text_tag.text:
            raw_text = text_tag.text
        else:
            raw_text = None

        # TODO: factor this out and reuse fix_agents
        db_refs = {}
        # Save raw text if available
        if raw_text:
            db_refs['TEXT'] = raw_text
        agent_name = raw_text
        # If we have a proper UID then we try to reconstruct an Agent from that
        if uid is not None and len(uid.split(':')) == 2:
            db_ns, db_id = uid.split(':')
            be_id = famplex_map.get((db_ns, db_id))
            if be_id:
                db_refs[db_ns] = db_id
                db_refs['FPLX'] = be_id
                agent_name = be_id
            elif db_ns in ['UP', 'Uniprot']:
                id_from_mnemonic = uniprot_client.get_id_from_mnemonic(db_id)
                if id_from_mnemonic:
                    db_id = id_from_mnemonic
                db_refs['UP'] = db_id
                hgnc_id = uniprot_client.get_hgnc_id(db_id)
                if hgnc_id:
                    db_refs['HGNC'] = hgnc_id
                    agent_name = hgnc_client.get_hgnc_name(hgnc_id)
                else:
                    gene_name = uniprot_client.get_gene_name(db_id)
                    if gene_name:
                        agent_name = gene_name
            elif db_ns == 'NCIT':
                db_refs['NCIT'] = db_id
                target = ncit_map.get(db_id)
                if target:
                    db_refs[target[0]] = target[1]
                    if target[0] == 'HGNC':
                        up_id = hgnc_client.get_uniprot_id(target[1])
                        agent_name = hgnc_client.get_hgnc_name(target[1])
                        if up_id:
                            db_refs['UP'] = up_id
                    elif target[0] == 'UP':
                        agent_name = uniprot_client.get_gene_name(target[1])
                        if agent_name:
                            hgnc_id = hgnc_client.get_hgnc_id(agent_name)
                            if hgnc_id:
                                db_refs['HGNC'] = hgnc_id
            elif db_ns == 'FA':
                db_refs['NXP'] = 'FA:' + db_id
            elif db_ns == 'XFAM':
                db_refs['PF'] = db_id.split('.')[0]
            elif db_ns == 'CHEBI':
                db_refs['CHEBI'] = 'CHEBI:' + db_id
            elif db_ns in ['GO', 'MESH', 'FPLX']:
                db_refs[db_ns] = db_id
            # Handle old BE mappings and add them as FPLX
            elif db_ns == 'BE':
                db_refs['FPLX'] = db_id
            elif db_ns in ['PR', 'CO', 'CVCL', 'EFO', 'ORPHANET']:
                db_refs[db_ns] = db_id
            else:
                logger.warning('Unknown database name space %s' % db_ns)
        if not agent_name:
            if raw_text is not None:
                agent_name = raw_text
            else:
                return None

        assert(agent_name)

        agent = Agent(agent_name, db_refs=db_refs)
        return agent

    @staticmethod
    def _get_amino_acid(ref):
        name_tag = ref.find("var/[@name='name']")
        if name_tag is not None:
            return name_tag.text
        return None

    @staticmethod
    def _get_site(site):
        residue = None
        position = None
        residue_tag = site.find("ref/[@category='residue-on-protein']")
        if residue_tag is not None:
            residue_aa_tag = residue_tag.find("var/[@name='amino-acid']/ref")
            if residue_aa_tag is not None:
                residue = SparserXMLProcessor._get_amino_acid(residue_aa_tag)
                position_tag = residue_tag.find("var/[@name='position']")
                if position_tag is not None:
                    position = position_tag.text.strip()
        if residue and residue.startswith('phospho'):
            residue = residue[7:]
            try:
                residue = get_valid_residue(residue)
            except InvalidResidueError:
                logger.error('Invalid residue found: %s' % residue)
                residue = None
        return residue, position

    def _get_evidence(self, text):
        ev = Evidence(source_api='sparser', pmid=self.pmid, text=text)
        return ev


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
