from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
from indra.databases import uniprot_client, hgnc_client
from indra.statements import Statement, Evidence

active_forms_files = ['sources/cure-active-forms-2017-01-30.txt',
                      'sources/r3-egfr-signaling-active-forms-2016-10-18.txt']

def read_jsons(fname):
    all_json = []
    with open(fname, 'rt') as fh:
        lin = fh.readlines()
    lin = [l.strip() for l in lin if l.strip()]
    for l in lin:
        lin_json = json.loads(l)
        all_json.append(lin_json[0])
    return all_json

def read_stmts(fname):
    jsons = read_jsons(fname)
    stmts = []
    for js in jsons:
        st = Statement.from_json(json.dumps(js))
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
        # Correct evidence
        if isinstance(st.evidence, Evidence):
            st.evidence = [st.evidence]
        for ev in st.evidence:
            ev.source_api = 'r3'
            ev.pmid = None
            ev.source_id = None
            ev.annotations = {}
            ev.epistemics = {}
        fix_protein_grounding(st.agent)
        # Correct supports/supported by
        st.supports = []
        st.supported_by = []
        st.belief = 1
    return stmts

def fix_protein_grounding(agent):
    if not agent.db_refs.get('TEXT'):
        agent.db_refs['TEXT'] = agent.name
    up_id = agent.db_refs.get('UP')
    if up_id:
        hgnc_symbol = uniprot_client.get_gene_name(up_id)
        hgnc_id = hgnc_client.get_hgnc_id(hgnc_symbol)
        if hgnc_id:
            agent.db_refs['HGNC'] = hgnc_id

if __name__ == '__main__':
    stmts = []
    for aff in active_forms_files:
        stmts = read_stmts(aff)
