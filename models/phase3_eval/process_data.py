from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import numpy
import pandas
from copy import copy
from random import shuffle
from collections import OrderedDict
from indra.databases import uniprot_client
from indra.databases import hgnc_client
from indra.literature import pubmed_client
from indra.statements import Agent, ModCondition

data_file = 'data/Korkut et al. Data 05122017.xlsx'
antibody_map_file = 'data/antibody_site_map.csv'
drug_grounding_file = 'data/drug_grounding.csv'

def read_data(fname=data_file):
    """Returns the data as a dictionary."""
    if fname.endswith('2017.xlsx'):
        skiprows1 = []
        skiprows2 = []
    else:
        skiprows1 = [0]
        skiprows2 = range(5)
    data = {}
    data['protein'] = pandas.read_excel(fname, sheetname='Protein Data',
                                        skiprows=skiprows1, index_col=None)
    data['phenotype'] = pandas.read_excel(fname,
                                          sheetname='Phenotype Data',
                                          skiprows=skiprows1, index_col=None)
    data['antibody'] = pandas.read_excel(fname,
                                         sheetname='Antibody Data',
                                         skiprows=skiprows2, index_col=None)
    data['prediction'] = pandas.read_excel(fname,
                                           sheetname='Prediction Targets',
                                           index_col=None)
    return data

def get_all_antibodies(data):
    ab_names = list(data['protein'].columns[2:])
    return ab_names

def get_phos_antibodies(data):
    ab_names = get_all_antibodies(data)
    ab_phos = []
    phospho_aa = ['S', 'T', 'Y']
    for abn in ab_names:
        for aa in phospho_aa:
            if abn.find('_p%s' % aa) != -1:
                ab_phos.append(abn)
                break
    return ab_phos

def get_antibody_genes(data, ab_name):
    """Return the UniProt IDs corresponding to a specific antibody."""
    df_filt = data['antibody'][data['antibody']['Protein Data ID'] == ab_name]
    up_ids = df_filt['UniProt ID'].values[0].split(',')
    return up_ids

def get_agent_from_upid(up_id):
    """Get an Agent based on a UniProt ID."""
    gene_name = uniprot_client.get_gene_name(up_id)
    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
    a = Agent(gene_name, db_refs={'HGNC': hgnc_id, 'UP': up_id})
    return a

def get_ras227_genes():
    ras227_file = '../../data/ras_pathway_proteins.csv'
    df = pandas.read_csv(ras227_file, sep='\t', index_col=None, header=None,
                         encoding='utf-8')
    gene_names = [x for x in df[0]]
    return gene_names

def get_all_gene_names(data, out_file='prior_genes.txt'):
    """Return all gene names corresponding to all ABs."""
    filt = pandas.notnull(data['antibody']['Protein Data ID'])
    data_filt = data['antibody'][filt]
    gene_names = data_filt['Gene Name']
    uniprot_ids = data_filt['UniProt ID']
    all_genes = set()
    invalid_genes = set()
    for gn, upid in zip(gene_names, uniprot_ids):
        # Some entries are lists of genes separated by commas
        # and we also strip off extra spaces
        names = [x.strip() for x in gn.split(',')]
        ids = [x.strip() for x in upid.split(',')]
        names_from_ids = [uniprot_client.get_gene_name(x) for x in ids]
        # Find invalid gene names
        for name in names:
            if not hgnc_client.get_hgnc_id(name):
                print('Invalid or deprecated gene symbol: %s' % name)
                invalid_genes.add(name)
        # Find inconsistent gene names and UniProt IDs
        if set(names) != set(names_from_ids):
            print('Inconsistent entries:')
            print('- Given gene names: %s' % ','.join(names))
            print('- Genes from uniprot IDs: %s' % ','.join(names_from_ids))
        # Add both the gene names and the gene names derived from UniProt IDs
        all_genes = all_genes.union(set(names)).union(set(names_from_ids))
    # Finally remove the invalid gene names
    all_genes = list(all_genes.difference(invalid_genes))
    # Add drug target genes
    drug_targets = get_drug_targets()
    for targets in drug_targets.values():
        all_genes += targets
    # Add other important genes, for now, the RAS pathway
    all_genes += get_ras227_genes()
    all_genes = sorted(list(set(all_genes)))
    print('%d genes in total' % len(all_genes))
    with open(out_file, 'wb') as fh:
        for gene in all_genes:
            fh.write(('%s\n' % gene).encode('utf-8'))
    return all_genes

def get_gene_pmids(genes, out_file='pmids.txt'):
    all_pmids = set()
    for gene in genes:
        print(gene)
        pmids = pubmed_client.get_ids_for_gene(gene)
        all_pmids = all_pmids.union(set(pmids))
    all_pmids = sorted(list(all_pmids))
    shuffle(all_pmids)
    with open(out_file, 'wb') as fh:
        for pmid in all_pmids:
            fh.write(('%s\n' % pmid).encode('utf-8'))
    return all_pmids

def get_drugs(data):
    drugs_col = \
        data['protein']['Sample Description (drug abbre. | dose or time-point)']
    drug_abbrevs = set()
    for cond in drugs_col:
        terms = cond.split(',')
        for term in terms:
            da, dose = term.split('|')
            drug_abbrevs.add(da)
    drug_abbrevs = sorted(list(drug_abbrevs))
    return drug_abbrevs

def get_drug_targets(fname=None):
    if not fname:
        fname = drug_grounding_file
    df = pandas.read_csv(fname, index_col=None, header=None, encoding='utf-8')
    abbrevs = df[1]
    target_upids = df[6]
    targets = {}
    for abb, tupid in zip(abbrevs, target_upids):
        targets[abb] = [uniprot_client.get_gene_name(ui)
                        for ui in tupid.split(',')]
    return targets

def get_single_drug_treatments(data):
    drugs_col = \
        data['protein']['Sample Description (drug abbre. | dose or time-point)']
    drug_tx = []
    for cond in drugs_col:
        terms = cond.split(',')
        if len(terms) == 1:
            drug_tx.append(cond)
    return drug_tx

def get_midas_data(data, pkn_abs, pkn_drugs, out_file='MD-korkut.csv'):
    drug_abbrevs = get_drugs(data)
    drug_abbrevs = [x for x in drug_abbrevs if x in pkn_drugs]
    print(drug_abbrevs)
    phospho_abs = get_phos_antibodies(data)
    phospho_abs = set(phospho_abs).intersection(set(pkn_abs))
    drug_cols = ['TR:%s:Drugs' % dr for dr in drug_abbrevs]
    print(drug_cols)
    all_values = []
    for row in data['protein'].iterrows():
        row = row[1]
        # Prepare row values
        values = OrderedDict()
        values['TR:SKMEL133:CellLine'] = 1
        for dc in drug_cols:
            values[dc] = 0
        # control TR columns will be 0 for all but CellLine
        values_control = copy(values)
        # Get drug conditions for row
        drug_cond = row[1]
        terms = drug_cond.split(',')
        for term in terms:
            da, dose = term.split('|')
            if da in drug_abbrevs:
                values['TR:%s:Drugs' % da] = dose
        # Get data conditions for row
        for ab_name in phospho_abs:
            values['DV:%s' % ab_name] = row[ab_name]
            values_control['DV:%s' % ab_name] = 1
            values['DA:%s' % ab_name] = 1440
            values_control['DA:%s' % ab_name] = 0
        all_values.append(values)
        all_values.append(values_control)
    df = pandas.DataFrame.from_records(all_values)
    # Rename columns in MIDAS the same way they are renamed in SIF
    def rename_df_columns(df):
        cols = df.columns.tolist()
        new_cols = []
        for c in cols:
            if c.find('AND') != -1:
                c = c.replace('AND', 'A_ND')
            if c.find('-') != -1:
                c = c.replace('-', '_')
            if c[3].isdigit():
                c = c[0:3] + 'abc_' + c[3:]
            new_cols.append(c)
        df.columns = new_cols
        return df
    df = rename_df_columns(df)
    df = df.drop_duplicates()
    df.to_csv(out_file, index=False)
    return df

def get_dynamic_range(data):
    ncol = len(data['protein'].columns)
    ranges = {}
    for i in range(2, ncol):
        vals = data['protein'].iloc[:,i]
        rangev = numpy.max(vals) - numpy.min(vals)
        ranges[data['protein'].columns[i]] = rangev
    return ranges

def get_average_dev(data):
    ncol = len(data['protein'].columns)
    devs = {}
    for i in range(2, ncol):
        vals = data['protein'].iloc[:,i]
        dev = numpy.mean(numpy.abs(1-vals))
        devs[data['protein'].columns[i]] = dev
    return devs

def get_max_dev(data):
    ncol = len(data['protein'].columns)
    devs = {}
    for i in range(2, ncol):
        vals = data['protein'].iloc[:,i]
        dev = numpy.mean(numpy.abs(1-vals))
        devs[data['protein'].columns[i]] = dev
    return devs

def get_antibody_map(data):
    phos_ab_map = get_phospho_antibody_map()
    ab_map = {}
    for _, row in data['antibody'].iterrows():
        ab_name = row['Protein Data ID']
        if ab_name in phos_ab_map:
            continue
        upids = row['UniProt ID'].split(',')
        for upid in upids:
            hgnc_symbol = uniprot_client.get_gene_name(upid)
            hgnc_id = hgnc_client.get_hgnc_id(hgnc_symbol)
            target = Agent(hgnc_symbol,
                           db_refs={'UP': upid,'HGNC': hgnc_id})
            try:
                ab_map[ab_name].append(target)
            except KeyError:
                ab_map[ab_name] = [target]
    ab_map.update(phos_ab_map)
    return ab_map


def get_phospho_antibody_map(fname=antibody_map_file):
    # First gather the annotations for the phosphosites
    df = pandas.read_csv(fname, index_col=None, sep=',', encoding='utf8')
    antibody_map = {}

    for _, row in df.iterrows():
        ps = row['phosphosite']
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
        mc = ModCondition('phosphorylation', residue, position)
        sub.mods = [mc]
        if ps in antibody_map:
            found = False
            for p in antibody_map[ps]:
                if p.name == sub.name and p.mods[0].residue == residue and \
                    p.mods[0].position == position:
                    found = True
                    break
            if not found:
                antibody_map[ps].append(sub)
        else:
            antibody_map[ps] = [sub]
    return antibody_map

if __name__ == '__main__':
    data = read_data(data_file)
    #gene_names = get_all_gene_names(data)
    #pmids = get_gene_pmids(gene_names)
