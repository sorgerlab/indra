from indra.util import _require_python3
import re
import copy
import json
import numpy
import pickle
import pandas
import itertools
import collections
from copy import deepcopy
from indra.statements import *
from indra.literature import pubmed_client
from indra.databases import hgnc_client, uniprot_client, cbio_client

rppa_file = 'data/TableS1-Split.xlsx'
rppa_pkl = 'data/TableS1-Split.pkl'
expression_file = 'data/Expression_Filtered.csv'
mutation_file = 'data/WES_variants_filtered.csv'
mutation_effect_file = 'data/mutation_effects.tsv'


def read_rppa_data(fname=rppa_pkl):
    """Return RPPA data as a dict of median/std DataFrames."""
    # If the filename passed is a pickle, just load it and return
    if fname.endswith('pkl'):
        print('Loading data from %s' % fname)
        with open(fname, 'rb') as fh:
            data = pickle.load(fh)
            return data

    # Otherwise just read the excel file
    data = {}
    for cell_line in cell_lines:
        print('Reading data for %s' % cell_line)
        data[cell_line] = {}
        # Read both the median and the std sheet for each cell line
        for data_type, postfix in (('median', ''), ('std', '-std')):
            # Handle unpredictable number of extra rows before the actual
            # header row.
            i = 0
            while True:
                df = pandas.read_excel(fname, sheetname=(cell_line + postfix),
                                       skiprows=i, header=0)
                if df.columns[0] == 'Drug':
                    break
                i += 1
            data[cell_line][data_type] = df
    return data


def read_variants(gene_names):
    """Return genetic variants reported by Sanger (provided locally)"""
    def read_aa_change(aa_change):
        match = re.match('p.([A-Z])(\d+)([A-Z])', aa_change)
        if not match:
            return None
        groups = match.groups()
        if len(groups) != 3:
            return None
        return groups

    def read_for_cell_line(cell_line):
        miss = {}
        nons = []
        df_filt = df[df.SAMPLE == cell_line]
        for _, row in df_filt.iterrows():
            # Check for valid/usable gene name
            if row['Gene'] not in gene_names:
                continue
            if row['Classification'] == 'missense':
                aa_change = read_aa_change(row['AA'])
                if not aa_change:
                    continue
                if row['Gene'] not in miss:
                    miss[row['Gene']] = [aa_change]
                else:
                    miss[row['Gene']].append(aa_change)
            elif row['Classification'] == 'nonsense':
                nons.append(row['Gene'])
        return miss, nons

    variants = {'missense': {}, 'nonsense': {}}
    df = pandas.read_csv(mutation_file)
    cell_line_map = {'SK-MEL-28': 'SKMEL28', 'WM-115': 'WM115',
            'MZ7-mel': 'MZ7MEL', 'MMAC-SF': 'MMACSF', 'RVH-421': 'RVH421',
            'C32': 'C32', 'LOXIMVI': 'LOXIMVI'}
    for cell_line_df, cell_line in cell_line_map.items():
        variants['missense'][cell_line], variants['nonsense'][cell_line] = \
            read_for_cell_line(cell_line_df)
    return variants


def read_ccle_variants(gene_names):
    """Return genetic variants reported by CCLE (via web service)"""
    cell_lines_db = [cl + '_SKIN' for cl in cell_lines]
    nons = cbio_client.get_ccle_mutations(gene_names, cell_lines_db, 'nonsense')
    miss = cbio_client.get_ccle_mutations(gene_names, cell_lines_db, 'missense')
    tmp = copy.deepcopy(nons)
    for cell_line, content in tmp.items():
        for gene, muts in content.items():
            if not muts:
                nons[cell_line].pop(gene, None)
    tmp = copy.deepcopy(miss)
    for cell_line, content in tmp.items():
        for gene, muts in content.items():
            if not muts:
                miss[cell_line].pop(gene, None)
                continue
            # Check for usable AA substitution
            grouped_muts = []
            for mut in muts:
                match = re.match('([A-Z])(\d+)([A-Z])', mut)
                if not match:
                    continue
                groups = match.groups()
                if len(groups) != 3:
                    continue
                grouped_muts.append(groups)
            miss[cell_line][gene] = grouped_muts

    variants = {}
    variants['missense'] = miss
    variants['nonsense'] = nons
    return variants


def read_mutation_effects():
    """Read mutation effects from PathwayCommons as ActiveForms."""
    df = pandas.read_csv(mutation_effect_file, sep='\t')
    stmts = []
    for _, row in df.iterrows():
        # Check for usable AA substitution
        match = re.match('([A-Z])(\d+)([A-Z])', row['Substitution'])
        if not match:
            continue
        groups = match.groups()
        if len(groups) != 3:
            continue
        aa_from, position, aa_to = groups
        # Make a grounded agent for the ActiveForm
        agent = agent_from_gene_name(row['Gene'])
        if not agent.db_refs:
            continue
        # Check if the reported effect is positive or negative
        if 'decreasing' in row['Effect']:
            is_active = False
        elif 'increasing' in row['Effect']:
            is_active = True
        else:
            continue
        # Make an ActiveForm Statement
        mc = MutCondition(position, aa_from, aa_to)
        agent.mutations = [mc]
        af = ActiveForm(agent, 'activity', is_active)
        stmts.append(af)
    return stmts


def analyze_cna_mrna(genes, cell_line):
    """Plot CNA vs mRNA reported by CCLE for given genes in a cell line."""
    import matplotlib.pyplot as plt
    cna = cbio_client.get_ccle_cna(genes, [cell_line])[cell_line]
    mrna = cbio_client.get_ccle_mrna(genes, [cell_line])[cell_line]
    plt.ion()
    plt.figure()
    for gene, cna_val in cna.items():
        if cna_val is None or gene not in mrna:
            continue
        mrna_val = mrna[gene]
        if mrna_val is None:
            continue
        plt.plot(cna_val, mrna_val, 'ro', alpha=0.5)
    plt.xlim([-2.1, 2.1])
    cna_vals = list(range(-2, 3))
    plt.title('%s cell line' % cell_line)
    plt.xticks(cna_vals, [str(c) for c in cna_vals])
    plt.ylabel('mRNA amounts reported by CCLE')
    plt.xlabel('DNA copy number alteration reported by CCLE')
    plt.show()


def _read_gene_list(path):
    df = pandas.read_csv(path, sep='\t', error_bad_lines=False, header=None)
    gene_list = df[0].tolist()
    return sorted(list(set(gene_list)))


def get_gene_names(data):
    """Return genes relevant to the data, RAS signaling, and the Fallahi paper.
    """
    gene_names = []
    antibody_list = []
    exclude_columns = ['Drug', 'Time (hr)', 'Concentration (uM)']
    for cl in cell_lines:
        antibody_list += [x for x in data[cl]['median']
                          if x not in exclude_columns]
    antibody_list = list(set(antibody_list))
    antibody_gene_list = []
    for ab in antibody_list:
        antibody_gene_list += list(antibody_map[ab].keys())
    ras_gene_list = _read_gene_list('../../data/ras_pathway_proteins.csv')
    msb2015_gene_list = \
            _read_gene_list('../../data/MFS_MSB_2015_gene_list.csv')
    msb2017_gene_list = \
            _read_gene_list('../../data/MFS_MSB_2017_gene_list.csv')
    gene_names = (ras_gene_list + msb2015_gene_list +
                  msb2017_gene_list + antibody_gene_list)
    gene_names = sorted(list(set(gene_names)))
    return gene_names


def get_gene_pmids(gene_names):
    """Return PMIDs for all genes of interest."""
    genes_pmid_list = []
    for gene in gene_names:
        genes_pmid_list += pubmed_client.get_ids_for_gene(gene)
    genes_pmid_list = list(set(genes_pmid_list))
    print('Found %d PMIDs for genes' % len(genes_pmid_list))
    return genes_pmid_list


def get_drug_pmids():
    """Return PMIDs for all the drugs and their synonyms."""
    drugs_pmid_list = []
    for drug_synonyms in drug_names.values():
        for drug_synonym in drug_synonyms:
            drugs_pmid_list += pubmed_client.get_ids(drug_synonym, retmax=5000)
    drugs_pmid_list = list(set(drugs_pmid_list))
    print('Found %d PMIDs for drugs' % len(drugs_pmid_list))
    return drugs_pmid_list


def get_all_pmids(data, pmid_file='pmids.txt'):
    """Return all PMIDs for genes and drugs of interest."""
    gene_names = get_gene_names(data)
    genes_pmid_list = get_gene_pmids(gene_names)
    drugs_pmid_list = get_drug_pmids()
    pmid_list = list(set(genes_pmid_list + drugs_pmid_list))
    with open(pmid_file, 'wb') as fh:
        for pmid in pmid_list:
            fh.write(('%s\n' % pmid).encode('utf-8'))
    return pmid_list


def find_extremes(data, fold_change, save_file=None):
    """Return rows of data which are above or below the given fold change."""
    liml, limu = (numpy.log2(1.0 / fold_change), numpy.log2(fold_change))
    all_extremes = []
    for cell_line in cell_lines:
        df = data[cell_line]['median']
        antibodies = df.columns[3:]
        for ab in antibodies:
            extremes = df.loc[(df[ab] < liml) | (df[ab] > limu)]
            for _, extreme in extremes.iterrows():
                drug = extreme['Drug']
                time = extreme['Time (hr)']
                conc = extreme['Concentration (uM)']
                val = extreme[ab]
                all_extremes.append([drug, time, conc, cell_line, ab, val])
    # Sort values
    all_extremes = sorted(all_extremes, key=lambda x: abs(x[5]), reverse=True)
    # Optionally save into a CSV file
    if save_file:
        with open(save_file, 'w') as fh:
            fh.write('Drug,Time (hr),Concentration (uM),' +
                     'CellLine,Antibody,Value\n')
            for vals in all_extremes:
                fh.write(','.join([str(v) for v in vals]) + '\n')
    return all_extremes


def find_cell_line_vars(data, fold_change, save_file=None):
    """Return conditions in which cell lines are qualitatively different."""
    liml, limu = (numpy.log2(1.0 / fold_change), numpy.log2(fold_change))
    all_vals = []
    # Look at all cell line combinations for differences
    for cl1, cl2 in itertools.combinations(cell_lines, 2):
        df1 = data[cl1]['median']
        df2 = data[cl2]['median']
        antibodies = df1.columns[3:]
        for ab in antibodies:
            filt = (((df1[ab] < liml) & (df2[ab] > limu)) |
                    ((df2[ab] < liml) & (df1[ab] > limu)))
            if filt.any():
                cl1_row = df1.loc[filt]
                cl2_row = df2.loc[filt]
                vals = [cl1_row['Drug'].values[0],
                        cl1_row['Time (hr)'].values[0],
                        cl1_row['Concentration (uM)'].values[0],
                        ab, cl1, cl2,
                        cl1_row[ab].values[0], cl2_row[ab].values[0]]
                all_vals.append(vals)
    # Sort values
    all_vals = sorted(all_vals, key=lambda x: abs(x[6] - x[7]), reverse=True)
    # Optionally save into a CSV file
    if save_file:
        with open(save_file, 'w') as fh:
            fh.write('Drug,Time (hr),Concentration (uM),Antibody,' +
                     'CellLine1,CellLine2,Value1,Value2\n')
            for vals in all_vals:
                fh.write(','.join([str(v) for v in vals]) + '\n')
    return all_vals


def agent_phos(name, phos_sites):
    """Return an INDRA agent from a name and list of phos sites."""
    agent = agent_from_gene_name(name)
    for residue, position in phos_sites:
        mc = ModCondition('phosphorylation', residue, position, True)
        agent.mods.append(mc)
    return agent


def agent_from_gene_name(name):
    """Return a grounded Agent based on a gene name."""
    agent = Agent(name)
    hgnc_id = hgnc_client.get_hgnc_id(name)
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    agent.db_refs = {'HGNC': hgnc_id, 'UP': uniprot_id}
    return agent


def get_antibody_agents():
    """Return a list of INDRA Agents corresponding to each antibody."""
    antibody_agents = {}
    for ab, agents in antibody_map.items():
        antibody_agents[ab] = []
        for name, phos_sites in agents.items():
            agent = agent_phos(name, phos_sites)
            antibody_agents[ab].append(agent)
    return antibody_agents


def get_agent_values_for_condition(data, cell_line, drug, time, dose):
    df = data[cell_line]['median']
    df = df[df['Drug'] == drug]
    df = df[df['Time (hr)'] == time]
    df = df[df['Concentration (uM)'] == dose]
    values = {}
    for ab in df.columns[3:]:
        value = df[ab].values[0]
        if not pandas.isnull(value):
            values[ab] = value
    return values


def build_context_json(data, cell_lines, filename="fallahi"):
    """Save context in a format readable by the pathway map"""
    select_df = pandas.DataFrame(columns=['Drug',
                                          'Concentration (uM)',
                                          'Time (hr)'])
    long_df = pandas.DataFrame(columns=['id', 'variable', 'value'])
    for cl in cell_lines:
        df = deepcopy(data[cl]['median'])
        df['cell_line'] = cl
        select_df = select_df.append(df[['Drug',
                                         'Concentration (uM)',
                                         'Time (hr)',
                                         'cell_line']])
        df['id'] = df.apply(lambda x: '%s_%s_%s_%s' %
                            (x['Drug'], x['Concentration (uM)'],
                             x['Time (hr)'], x['cell_line']),
                            axis=1)
        df = df.drop(['Drug', 'Time (hr)',
                      'Concentration (uM)', 'cell_line'], axis=1)
        df = df.melt(id_vars=['id'])
        long_df = long_df.append(df)
    ab_dict = {'antibody': [], 'gene': [], 'site': []}
    for ab, d in antibody_map.items():
        for hgnc, sites in d.items():
            if len(sites) > 0:
                for s in sites:
                    residue = s[0] + s[1]
                    ab_dict['antibody'].append(ab)
                    ab_dict['gene'].append(hgnc)
                    ab_dict['site'].append(residue)
            else:
                residue = None
                ab_dict['antibody'].append(ab)
                ab_dict['gene'].append(hgnc)
                ab_dict['site'].append(residue)
    ab_df = pandas.DataFrame(ab_dict)
    long_data = pandas.merge(long_df, ab_df,
                             left_on='variable', right_on='antibody')
    long_data = long_data[~long_data['value'].isin([None])]
    recursivedict = lambda: collections.defaultdict(recursivedict)
    j = recursivedict()
    df1 = long_data
    ids = df1['id'].unique().tolist()
    for i in ids:
        df2 = df1[df1['id'] == i]
        genes = df2['gene'].unique().tolist()
        for g in genes:
            df3 = df2[df2['gene'] == g]
            df4 = df3[~df3['site'].isin([None])]
            j[i][g]['members']['antibodies'] = df4['antibody'].tolist()
            j[i][g]['members']['sites'] = df4['site'].tolist()
            j[i][g]['members']['values'] = df4['value'].tolist()
            df5 = df3[df3['site'].isin([None])]
            j[i][g]['node']['antibodies'] = df5['antibody'].tolist()
            j[i][g]['node']['values'] = df5['value'].tolist()
    with open(filename + '_data.json', 'w') as outfile:
        json.dump(j, outfile, indent=None, sort_keys=True)
    select_df = select_df.drop_duplicates()
    select_df = select_df[['Drug', 'Concentration (uM)',
                           'Time (hr)', 'cell_line']]
    select_array = select_df.values.tolist()
    with open(filename + '_select.json', 'w') as outfile:
        json.dump(select_array, outfile, indent=None, sort_keys=False)


cell_lines = ['C32', 'LOXIMVI', 'RVH421', 'SKMEL28', 'WM115',
              'COLO858', 'K2', 'MMACSF', 'MZ7MEL', 'WM1552C']

drug_names = {'AZ628': ['AZ628', 'AZ_628', 'AZ-628', '878739-06-1'],
             'Selumetinib': ['Selumetinib', 'AZD6244', 'AZD 6244',
                             'ARRY-142886', '606143-52-6'],
             'SB590885': ['SB590885', 'SB-590885', 'J-501805', '405554-55-4'],
             'Vemurafenib': ['Vemurafenib', 'Zelboraf', 'PLX4032', 'PLX-4032',
                             'RG7204', 'RG-7204', 'R05185426', '918504-65-1',
                             '1029872-54-5'],
             'PLX4720': ['PLX4720', 'PLX 4720', 'PLX-4720', 'PLX_4720']}


drug_targets = {'AZ628': ['BRAF'],
                'Selumetinib': ['MAP2K1', 'MAP2K2'],
                'SB590885': ['BRAF'],
                'Vemurafenib': ['BRAF'],
                'PLX4720': ['BRAF']}

drug_grounding = {
        'AZ628': {'CHEBI': '1'},
        'Selumetinib': {'CHEBI': '2'},
        'SB590885': {'CHEBI': '3'},
        'Vemurafenib': {'CHEBI': '4'},
        'PLX4720': {'CHEBI': '5'},
        }

drug_doses = [0.0032, 0.01, 0.0316, 0.1, 0.316, 1.0, 3.16]



antibody_map = {
    'pMEK(S217/221)':
        {'MAP2K1': [('S', '218')],
         'MAP2K2': [('S', '222')]},
    'pERK(T202/Y204)':
        {'MAPK1': [('T', '185'), ('Y', '187')],
         'MAPK3': [('T', '202'), ('Y', '204')]},
    'p-p90RSK(S380)':
        {'RPS6KA1': [('S', '380')]},
    'p-p90RSK(T573)':
        {'RPS6KA1': [('T', '573')]},
    'p-AKT(T308)':
        {'AKT1': [('T', '308')]},
    'p-AKT(S473)':
        {'AKT1': [('S', '473')],
         'AKT2': [('S', '474')],
         'AKT3': [('S', '472')]},
    'p-mTOR(S2448)':
        {'MTOR': [('S', '2448')]},
    'p-p70S6K(T421/S424)':
        {'RPS6KB1': [('T', '444'), ('S', '447')]},
    'p-p70S6K(T389)':
        {'RPS6KB1': [('T', '390')]},
    'p-S6(S235/236)':
        {'RPS6': [('S', '235'), ('S', '236')]},
    'p-AMPK(T172)':
        {'PRKAA1': [('T', '183')]},
    'p-JNK(T183/Y185)':
        {'MAPK8': [('T', '183'), ('Y', '185')],
         'MAPK9': [('T', '183'), ('Y', '185')],
         'MAPK10': [('T', '221'), ('Y', '223')]},
    'Total c-Jun':
        {'JUN': []},
    'p-cJun(S63)':
        {'JUN': [('S', '63')]},
    'p-P38(T180/Y182)':
        {'MAPK14': [('T', '180'), ('Y', '182')]},
    'p-HSP27(S82)':
        {'HSPB1': [('S', '82')]},
    'p-NFKB(S536)':
        {'RELA': [('S', '536')]},
    'Bim':
        {'BCL2L11': []},
    'cPARP':
        {'PARP1': []},
    'p-Histone H3(S10)':
        {'HIST1H3A': [('S', '11')],
         'HIST1H3B': [('S', '11')],
         'HIST1H3C': [('S', '11')],
         'HIST1H3D': [('S', '11')],
         'HIST1H3E': [('S', '11')],
         'HIST1H3F': [('S', '11')],
         'HIST1H3G': [('S', '11')],
         'HIST1H3H': [('S', '11')],
         'HIST1H3I': [('S', '11')],
         'HIST1H3J': [('S', '11')],
         'HIST2H3A': [('S', '11')],
         'HIST2H3C': [('S', '11')],
         'HIST2H3D': [('S', '11')],
         'H3F3A': [('S', '11')],
         'H3F3B': [('S', '11')]},
    'p27 Kip1':
        {'CDKN1B': []}
    }
