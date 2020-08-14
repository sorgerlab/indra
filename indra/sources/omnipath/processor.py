from __future__ import unicode_literals
import logging
from collections import Counter
from indra.preassembler.grounding_mapper import GroundingMapper
from indra.statements import modtype_to_modclass, Agent, Evidence, Complex

logger = logging.getLogger(__file__)


class OmniPathProcessor(object):
    def __init__(self, ptm_json=None, ligrec_json=None):
        self.statements = []
        self.statements.extend(self._stmts_from_op_mods(ptm_json))
        self.statements.extend(self._stmt_from_pp_lr(ligrec_json))

    def _stmts_from_op_mods(self, ptm_json):
        """Build Modification Statements from a list of Omnipath PTM entries
        """
        ptm_stmts = []
        unhandled_mod_types = []
        if ptm_json is None:
            return []
        for mod_entry in ptm_json:
            enz = self._agent_from_up_id(mod_entry['enzyme'])
            sub = self._agent_from_up_id(mod_entry['substrate'])
            res = mod_entry['residue_type']
            pos = mod_entry['residue_offset']
            # evidence = [Evidence('omnipath', None, pmid)
            #            for pmid in mod_entry['references']]
            evidence = [Evidence('omnipath', None, None)]
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

    def _stmt_from_pp_lr(self, ligrec_json):
        """Build Complex statements from a json containing ligand-receptor
        interactions built from PyPath
        """
        stmt_list = []
        if ligrec_json is None:
            return stmt_list

        for lr_entry in ligrec_json:
            for source_pmid in lr_entry['references']:
                source_id, pmid = source_pmid.split(':')
                agents = self._complex_agents_from_op_complex(
                    [lr_entry['source'], lr_entry['target']]
                )
                evidence = Evidence(source_api='omnipath',
                                    source_id=source_id,
                                    pmid=pmid,
                                    annotations=None,
                                    text_refs=None)
                stmt_list.append(Complex(members=agents,
                                         evidence=evidence))

        return stmt_list

    @staticmethod
    def _agent_from_up_id(up_id):
        """Build an Agent object from a Uniprot ID. Adds db_refs for both
        Uniprot and HGNC where available."""
        db_refs = {'UP': up_id}
        ag = Agent(None, db_refs=db_refs)
        GroundingMapper.standardize_agent_name(ag)
        return ag

    def _complex_agents_from_op_complex(self, up_id_list):
        """Return a list of agents from a string containing multiple UP ids
        """
        # Return list of contained agents
        if isinstance(up_id_list, list):
            return [self._agent_from_up_id(up_id) for up_id in up_id_list]
        elif isinstance(up_id_list, str):
            return [self._agent_from_up_id(up_id_list)]
        else:
            raise TypeError('Unable to produce agents from object %s' %
                            up_id_list.__class__)
