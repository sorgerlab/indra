from indra.util import _require_python3
import numpy
import pandas
import itertools
from indra.literature import pubmed_client
from indra.statements import Agent, ModCondition
from indra.databases import hgnc_client, uniprot_client

rppa_file = 'data/TableS1-Split.xlsx'
expression_file = 'data/Expression_Filtered.csv'
mutation_file = 'data/WES_variants_filtered.csv'


def read_rppa_data(fname=rppa_file):
    """Return RPPA data as a dict median/std DataFrames."""
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
            _read_gene_list('../../data/MohammadFS_MSB_2015_gene_list.csv')
    msb2017_gene_list = \
            _read_gene_list('../../data/MohammadFS_MSB_2017_gene_list.csv')
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
    for drug_synonyms in drug_dict.values():
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
    agent = Agent(name)
    hgnc_id = hgnc_client.get_hgnc_id(name)
    uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
    agent.db_refs = {'HGNC': hgnc_id, 'UP': uniprot_id}
    for residue, position in phos_sites:
        mc = ModCondition('phosphorylation', residue, position, True)
        agent.mods.append(mc)
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


cell_lines = ['C32', 'COLO858', 'K2', 'LOXIMVI', 'MMACSF', 'MZ7MEL',
              'RVH421', 'SKMEL28', 'WM115', 'WM1552C']

drug_dict = {'AZ628': ['AZ628', 'AZ_628', 'AZ-628', '878739-06-1'],
             'Selumetinib': ['Selumetinib', 'AZD6244', 'AZD 6244',
                             'ARRY-142886', '606143-52-6'],
             'SB590885': ['SB590885', 'SB-590885', 'J-501805', '405554-55-4'],
             'Vemurafenib': ['Vemurafenib', 'Zelboraf', 'PLX4032', 'PLX-4032',
                             'RG7204', 'RG-7204', 'R05185426', '918504-65-1',
                             '1029872-54-5'],
             'PLX4720': ['PLX4720', 'PLX 4720', 'PLX-4720', 'PLX_4720']}

antibody_map = {
    'pMEK(S217/221)':
        {'MAP2K1': [('S', '218')],
         'MAP2K2': [('S', '222')]},
    'pERK(T202/Y204)':
        {'MAPK1': [('T', '202'), ('Y', '204')],
         'MAPK3': [('T', '185'), ('Y', '187')]},
    'p-p90RSK(S380)':
        {'RPS6KA1': [('S', '380')]},
    'p-p90RSK(T573)':
        {'RPS6KA1': [('T', '573')]},
    'p-AKT(T308)':
        {'AKT1': [('T', '308')]},
    'p-AKT(S473)':
        {'AKT1': [('S', '473')],
         'AKT2': [('S', '474')],
         'AKT3': [('S', '475')]},
    'p-mTOR(S2448)':
        {'MTOR': [('S', '2448')]},
    'p-p70S6K(T421/S424)':
        {'RPS6KB1': [('T', '421'), ('S', '424')]},
    'p-p70S6K(T389)':
        {'RPS6KB1': [('T', '389')]},
    'p-S6(S235/236)': 
        {'RPS6': [('S', '235'), ('S', '236')]},
    'p-AMPK(T172)':
        {'PRKAA1': [('T', '172')]},
    'p-JNK(T183/Y185)':
        {'MAPK8': [('T', '183'), ('Y', '185')],
         'MAPK9': [('T', '183'), ('Y', '185')],
         'MAPK10': [('T', '183'), ('Y', '185')]},
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
        {'HIST1H3A': [('S', '10')],
         'HIST1H3B': [('S', '10')],
         'HIST1H3C': [('S', '10')],
         'HIST1H3D': [('S', '10')],
         'HIST1H3E': [('S', '10')],
         'HIST1H3F': [('S', '10')],
         'HIST1H3G': [('S', '10')],
         'HIST1H3H': [('S', '10')],
         'HIST1H3I': [('S', '10')],
         'HIST1H3J': [('S', '10')],
         'HIST2H3A': [('S', '10')],
         'HIST2H3C': [('S', '10')],
         'HIST2H3D': [('S', '10')],
         'H3F3A': [('S', '10')],
         'H3F3B': [('S', '10')]},
    'p27 Kip1':
        {'CDKN1B': []}
    }



if __name__ == '__main__':
    rppa_data = read_rppa_data(rppa_file)
    all_pmids = get_all_pmids(rppa_data, pmid_file='pmids.txt')
