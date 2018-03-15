from indra.util import _require_python3
import json
import pickle
from indra.databases import uniprot_client, hgnc_client
from indra.statements import stmts_from_json, Evidence
from util import prefixed_pkl

active_forms_files = ['cure-active-forms-2017-03-30.txt']


def read_jsons(fname):
    all_json = []
    with open(fname, 'rt') as fh:
        lin = fh.readlines()
    lin = [l.strip() for l in lin if l.strip()]
    for l in lin:
        lin_json = json.loads(l)
        all_json.append(lin_json)
    return all_json


def read_stmts(fname):
    jsons = read_jsons(fname)
    stmts = []
    for js in jsons:
        if not isinstance(js['evidence'], list):
            js['evidence'] = [js['evidence']]
        st = stmts_from_json([js])[0]
        if not hasattr(st.agent, 'mutations'):
            st.agent.mutations = []
        if not hasattr(st.agent, 'location'):
            st.agent.location = None
        if not hasattr(st.agent, 'activity'):
            st.agent.activity = None
        # Location and bound conditions alone not relevant here
        if (not st.agent.mods) and (not st.agent.mutations):
            continue
        st.agent.location = None
        st.agent.bound_conditions = []
        stmts.append(st)
    for st in stmts:
        for ev in st.evidence:
            ev.pmid = None
            ev.source_id = None
            ev.annotations = {}
            ev.epistemics = {}
        fix_mod_conditions(st.agent)
        fix_protein_grounding(st.agent)
        # Correct supports/supported by
        st.supports = []
        st.supported_by = []
        st.belief = 1
    stmts = eliminate_duplicates(stmts)
    return stmts


def fix_mod_conditions(agent):
    for mc in agent.mods:
        if mc.mod_type == 'glucosylation':
            mc.mod_type = 'glycosylation'
        if mc.mod_type == 'fucosylation':
            mc.mod_type = 'glycosylation'


def fix_protein_grounding(agent):
    for k, v in agent.db_refs.items():
        agent.db_refs.pop(k, None)
        agent.db_refs[k.upper()] = v
    if not agent.db_refs.get('TEXT'):
        agent.db_refs['TEXT'] = agent.name
    up_id = agent.db_refs.get('UP')
    if up_id:
        up_id = up_id.split('-')[0]
        agent.db_refs['UP'] = up_id
        hgnc_symbol = uniprot_client.get_gene_name(up_id)
        hgnc_id = hgnc_client.get_hgnc_id(hgnc_symbol)
        if hgnc_id:
            agent.name = hgnc_symbol
            agent.db_refs['HGNC'] = hgnc_id


def eliminate_duplicates(stmts):
    stmts_by_ag = {}
    for stmt in stmts:
        stmt.agent.db_refs.pop('TEXT', None)
        try:
            stmts_by_ag[stmt.agent.name].append(stmt)
        except KeyError:
            stmts_by_ag[stmt.agent.name] = [stmt]
    unique_by_ag = {}
    for k, v in stmts_by_ag.items():
        for st in v:
            found = False
            try:
                uniques = unique_by_ag[k]
            except KeyError:
                unique_by_ag[k] = []
                uniques = []
            for stmt in uniques:
                if stmt.equals(st):
                    found = True
                    break
            if not found:
                unique_by_ag[k].append(st)
    new_stmts = []
    for k, v in unique_by_ag.items():
        new_stmts += v
    return new_stmts

if __name__ == '__main__':
    stmts = []
    for aff in active_forms_files:
        stmts = read_stmts(aff)
    with open(prefixed_pkl('r3'), 'wb') as fh:
        pickle.dump(stmts, fh)
