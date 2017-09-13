import numpy
import pandas
import itertools

rppa_file = 'data/TableS1-Split.xlsx'
expression_file = 'data/Expression_Filtered.csv'
mutation_file = 'data/WES_variants_filtered.csv'

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

antibody_HGNC_dict = {'pMEK(S217/221)': ['MAP2K1', 'MAP2K2'],
                      'pERK(T202/Y204)': ['MAPK1', 'MAPK3'],
                      'p-p90RSK(S380)': ['RPS6KA1'],
                      'p-p90RSK(T573)': ['RPS6KA1'],
                      'p-AKT(T308)': ['AKT1'],
                      'p-AKT(S473)': ['AKT1', 'AKT2', 'AKT3'],
                      'p-mTOR(S2448)': ['MTOR'],
                      'p-p70S6K(T421/S424)': ['RPS6KB1'],
                      'p-p70S6K(T389)': ['RPS6KB1'],
                      'p-S6(S235/236)': ['RPS6'],
                      'p-AMPK(T172)': ['PRKAA1'],
                      'p-JNK(T183/Y185)': ['MAPK8', 'MAPK9', 'MAPK10'],
                      'Total c-Jun': ['JUN'],
                      'p-cJun(S63)': ['JUN'],
                      'p-P38(T180/Y182)': ['MAPK14'],
                      'p-HSP27(S82)': ['HSPB1'],
                      'p-NFKB(S536)': ['RELA'],
                      'Bim': ['BCL2L11'],
                      'cPARP': ['PARP1'],
                      'p-Histone H3(S10)': ['HIST1H3A', 'HIST1H3B', 'HIST1H3C',
                                            'HIST1H3D', 'HIST1H3E', 'HIST1H3F',
                                            'HIST1H3G', 'HIST1H3H', 'HIST1H3I',
                                            'HIST1H3J', 'HIST2H3A', 'HIST2H3C',
                                            'HIST2H3D', 'H3F3A', 'H3F3B'],
                      'p27 Kip1': ['CDKN1B']}


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


def get_gene_names(data):
    # Get the gene names in the data from the data frame
    # Append any other gene names that are also relevant that are not in
    #   the data
    # Return gene_names
    return None


def get_gene_pmids(gene_names):
    # Use indra.literature.pubmed_client.get_ids_for_gene for each gene
    #   to get its PMIDs
    # Collect all PMIDs in one list and return
    return None


def get_drug_pmids(data):
    # Extract name of drugs from the data or make a list at the top of
    #   this file
    # Use indra.literature.pubmed_client.get_ids with the
    #   name and/or synonym(s) of each drug and get the relevant PMIDS
    # Collect all PMIDs in one list and return
    return None


def get_all_pmids(data, pmid_file='pmids.txt'):
    # Call both get_gene_pmids and get_drug_pmids, and combine the results
    # Save list of PMIDs in the pmid_file path given
    # Return all PMIDs in a single list
    return None


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
