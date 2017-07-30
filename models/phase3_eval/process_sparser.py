from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import glob
import json
from indra.sources import sparser
from indra.statements import *
from indra.databases import uniprot_client, hgnc_client

base_folder = 'sources/sparser-20170530'

def get_file_names(base_dir):
    fnames = glob.glob(os.path.join(base_dir, '*.json'))
    return fnames

def get_file_stmts(fname):
    with open(fname, 'rt') as fh:
        print(fname)
        try:
            jd = json.load(fh)
        except ValueError as e:
            print(e)
            return []
    for st in jd:
        if st.get('type') == 'Translocation':
            for loc in ['from_location', 'to_location']:
                val = st.get(loc)
                try:
                    loc_valid = get_valid_location(val)
                    st[loc] = loc_valid
                except:
                    st[loc] = None
        try:
            res = st['residue']
            if res is False:
                st['residue'] = None
        except:
            pass

        try:
            res = st.get('residue')
            if res:
                get_valid_residue(res)
        except:
            st['residue'] = None

        try:
            res = st['position']
            if res is False:
                st['position'] = None
        except:
            pass

    stmts = stmts_from_json(jd)
    stmts = fix_stmts(stmts)
    return stmts

def read_stmts(folder):
    fnames = get_file_names(folder)
    all_stmts = []
    for fname in fnames:
        st = get_file_stmts(fname)
        all_stmts += st
    return all_stmts

def fix_stmts(stmts):
    new_stmts = []
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid and ev.pmid.startswith('PMID'):
                ev.pmid = ev.pmid[:-4]
        # Skip if no subject
        if isinstance(stmt, RegulateActivity):
            if stmt.subj is None:
                continue
        # Skip if no locations
        if isinstance(stmt, Translocation):
            if not (stmt.from_location or stmt.to_location):
                continue
        for agent in stmt.agent_list():
            if agent is not None:
                upid = agent.db_refs.get('UP')
                if upid:
                    gene_name = uniprot_client.get_gene_name(upid)
                    if gene_name:
                        agent.name = gene_name
                        if uniprot_client.is_human(upid):
                            hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                            if hgnc_id:
                                agent.db_refs['HGNC'] = hgnc_id

        new_stmts.append(stmt)
    return new_stmts

if __name__ == '__main__':
    stmts = read_stmts(base_folder)
