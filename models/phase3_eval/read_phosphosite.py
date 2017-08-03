from __future__ import print_function, absolute_import, unicode_literals
from builtins import dict, str
import pandas
from copy import deepcopy
from indra.sources import biopax
from indra.util import unicode_strs
from indra.sources import biopax
from indra.databases import hgnc_client, uniprot_client
from indra.statements import Agent, Phosphorylation, Evidence

phosphosite_file = 'sources/annotated_kinases_v5.csv'
phosphosite_owl_file = 'sources/Kinase_substrates.owl'

def read_phosphosite_owl(fname=phosphosite_owl_file):
    bp = biopax.process_owl(fname)
    for stmt in bp.statements:
        for ev in stmt.evidence:
            ev.source_api = 'phosphosite'
            ev.epistemics = {'direct': True}
    return bp.statements

def read_phosphosite(fname=phosphosite_file):
    df = pandas.read_csv(fname, index_col=None, sep='\t', encoding='utf8')
    statements = []
    for _, row in df.iterrows():
        sub_upid = row['SUB_ID']
        if not pandas.isnull(sub_upid):
            if sub_upid.find('-') != -1:
                sub_upid = sub_upid.split('-')[0]
            sub_hgnc_symbol = uniprot_client.get_gene_name(sub_upid)
            sub_hgnc = hgnc_client.get_hgnc_id(sub_hgnc_symbol)
        else:
            sub_hgnc_symbol = row['SUB_GENE']
            sub_hgnc_id = hgnc_client.get_hgnc_id(sub_hgnc_symbol)
            sub_upid = hgnc_client.get_uniprot_id(sub_hgnc_id)
            if sub_upid is None:
                continue
        sub = Agent(sub_hgnc_symbol,
                    db_refs={'UP': sub_upid,'HGNC': sub_hgnc})
        residue = row['Actual_site'][0]
        if len(row['Actual_site']) > 1:
            position = row['Actual_site'][1:]
        else:
            position = None

        ps = row['phosphosite']
        kin_upid = row['KIN_ID']
        if not pandas.isnull(kin_upid):
            if kin_upid.find('-') != -1:
                kin_upid = kin_upid.split('-')[0]
            if not uniprot_client.is_human(kin_upid):
                continue
            kin_hgnc_symbol = uniprot_client.get_gene_name(kin_upid)
            kin_hgnc = hgnc_client.get_hgnc_id(kin_hgnc_symbol)
        else:
            kin_hgnc_symbol = row['KINASE_GENE_SYMBOL']
            if pandas.isnull(kin_hgnc_symbol):
                continue
            kin_hgnc_id = hgnc_client.get_hgnc_id(kin_hgnc_symbol)
            kin_upid = hgnc_client.get_uniprot_id(kin_hgnc_id)
            if kin_upid is None:
                continue
        kin = Agent(kin_hgnc_symbol,
                    db_refs={'UP': kin_upid, 'HGNC': kin_hgnc})

        ev = Evidence(source_api='phosphosite')
        st = Phosphorylation(kin, sub, residue, position, evidence = [ev])
        statements.append(st)
    return statements
