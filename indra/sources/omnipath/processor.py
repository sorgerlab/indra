from __future__ import unicode_literals
import logging
from collections import Counter
from indra.ontology.standardize import standardize_agent_name
from indra.statements import modtype_to_modclass, Agent, Evidence, Complex, \
    get_statement_by_name as stmt_by_name, BoundCondition

logger = logging.getLogger(__name__)


class OmniPathProcessor(object):
    def __init__(self, ptm_json=None, ligrec_json=None):
        self.statements = []
        self.ptm_json = ptm_json
        self.ligrec_json = ligrec_json

    def process_ptm_mods(self):
        """Process ptm json if present"""
        if self.ptm_json:
            self.statements += self._stmts_from_op_mods(self.ptm_json)

    def process_ligrec_interactions(self):
        """Process ligand-receptor json if present"""
        if self.ligrec_json:
            self.statements += self._stmt_from_op_lr(self.ligrec_json)

    def _stmts_from_op_mods(self, ptm_json):
        """Build Modification Statements from a list of Omnipath PTM entries
        """
        ptm_stmts = []
        unhandled_mod_types = []
        annot_ignore = {'enzyme', 'substrate', 'residue_type',
                        'residue_offset', 'references', 'modification'}
        if ptm_json is None:
            return []
        for mod_entry in ptm_json:
            enz = self._agent_from_up_id(mod_entry['enzyme'])
            sub = self._agent_from_up_id(mod_entry['substrate'])
            res = mod_entry['residue_type']
            pos = mod_entry['residue_offset']
            evidence = []
            for source_pmid in mod_entry['references']:
                source_db, pmid = source_pmid.split(':', 1)
                if 'pmc' in pmid.lower():
                    text_refs = {'PMCID': pmid.split('/')[-1]}
                    pmid = None
                else:
                    text_refs = None
                evidence.append(Evidence(
                    source_api='omnipath',
                    source_id=source_db,
                    pmid=pmid,
                    text_refs=text_refs,
                    annotations={k: v for k, v in mod_entry.items() if k not
                                 in annot_ignore}
                ))
            mod_type = mod_entry['modification']
            modclass = modtype_to_modclass.get(mod_type)
            if modclass is None:
                unhandled_mod_types.append(mod_type)
                continue
            else:
                stmt = modclass(enz, sub, res, pos, evidence)
            ptm_stmts.append(stmt)
        print(Counter(unhandled_mod_types))
        return ptm_stmts

    def _stmt_from_op_lr(self, ligrec_json):
        """Make ligand-receptor Complexes from Omnipath API interactions db"""
        ligrec_stmts = []
        ign_annot = {'source_sub_id', 'source', 'target', 'references'}
        no_refs = 0
        bad_pmid = 0
        no_consensus = 0
        if ligrec_json is None:
            return ligrec_stmts

        for lr_entry in ligrec_json:
            if not lr_entry['references']:
                no_refs += 1
                continue
            if len(lr_entry['sources']) == 1 and \
                    lr_entry['sources'][0].lower() == 'protmapper':
                logger.warning('Protmapper only source, skipping...')
                continue

            # Assemble evidence
            evidence = []
            for source_pmid in lr_entry['references']:
                source_db, pmid = source_pmid.split(':')
                if len(pmid) > 8:
                    bad_pmid += 1
                    continue
                annot = {k: v for k, v in lr_entry.items() if k not in
                         ign_annot}
                annot['source_sub_id'] = source_db
                evidence.append(Evidence(source_api='omnipath', pmid=pmid,
                                         annotations=annot))

            # Get complexes
            ligrec_stmts.append(self._get_op_complex(lr_entry['source'],
                                                     lr_entry['target'],
                                                     evidence))

            # On consensus, make Activations or Inhibitions as well
            if bool(lr_entry['consensus_stimulation']) ^ \
               bool(lr_entry['consensus_inhibition']):
                activation = True if lr_entry['consensus_stimulation'] else \
                    False
                ligrec_stmts.append(self._get_ligrec_regs(
                    lr_entry['source'], lr_entry['target'], evidence,
                    activation=activation))
            elif lr_entry['consensus_stimulation'] and \
                    lr_entry['consensus_inhibition']:
                no_consensus += 1

        if no_refs:
            logger.warning(f'{no_refs} entries without references were '
                           f'skipped')
        if bad_pmid:
            logger.warning(f'{bad_pmid} references with bad pmids were '
                           f'skipped')
        if no_consensus:
            logger.warning(f'{no_consensus} entries with conflicting '
                           f'regulation were skipped')

        return ligrec_stmts

    @staticmethod
    def _agent_from_up_id(up_id):
        """Build an Agent object from a Uniprot ID. Adds db_refs for both
        Uniprot and HGNC where available."""
        db_refs = {'UP': up_id}
        ag = Agent(up_id, db_refs=db_refs)
        standardize_agent_name(ag)
        return ag

    def _bc_agent_from_up_list(self, up_id_list):
        # Return the first agent with the remaining agents as a bound condition
        agents_list = [self._agent_from_up_id(up_id) for up_id in up_id_list]
        agent = agents_list[0]
        agent.bound_conditions = \
            [BoundCondition(a, True) for a in agents_list[1:]]
        return agent

    def _complex_agents_from_op_complex(self, up_id_str):
        """Return a list of agents from a string containing multiple UP ids
        """
        # Get agents
        if 'complex' in up_id_str.lower():
            up_id_list = [up for up in up_id_str.split(':')[1].split('_')]
        else:
            up_id_list = [up_id_str]

        return [self._agent_from_up_id(up_id) for up_id in up_id_list]

    def _get_op_complex(self, source, target, evidence_list):
        ag_list = self._complex_agents_from_op_complex(source) + \
                  self._complex_agents_from_op_complex(target)
        return Complex(members=ag_list,
                       evidence=evidence_list)

    def _get_ligrec_regs(self, source, target, evidence_list, activation=True):
        # Check if any of the agents is a complex
        # Source
        if 'complex' in source.lower():
            # Make bound condition agent
            up_id_list = [up for up in source.split(':')[1].split('_')]
            subj = self._bc_agent_from_up_list(up_id_list)
        else:
            subj = self._agent_from_up_id(source)
        # Target
        if 'complex' in target.lower():
            # Make bound condition agent
            up_id_list = [up for up in target.split(':')[1].split('_')]
            obj = self._bc_agent_from_up_list(up_id_list)
        else:
            obj = self._agent_from_up_id(target)

        # Regular case:
        Regulation = stmt_by_name('activation') if activation else \
            stmt_by_name('inhibition')

        regulation = Regulation(subj=subj, obj=obj, evidence=evidence_list)
        return regulation
