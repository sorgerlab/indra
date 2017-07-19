import numpy
import pandas
import itertools

rppa_file = 'data/TableS1-Split.xlsx'
expression_file = 'data/Expression_Filtered.csv'
mutation_file = 'data/WES_variants_filtered.csv'

cell_lines = ['C32', 'COLO858', 'K2', 'LOXIMVI', 'MMACSF', 'MZ7MEL',
              'RVH421', 'SKMEL28', 'WM115', 'WM1552C']

def read_rppa_data(fname=rppa_file):
    """Return RPPA data as a dict median/std DataFrames."""
    data = {}
    for cell_line in cell_lines:
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


def find_extremes(data, fold_change):
    """Return rows of data which are above or below the given fold change."""
    liml, limu = (numpy.log2(1.0/fold_change), numpy.log2(fold_change))
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
                all_extremes.append([cell_line, ab, drug, time, conc, val])
    return all_extremes

def find_cell_line_vars(data, fold_change):
    """Return conditions in which cell lines are qualitatively different."""
    liml, limu = (numpy.log2(1.0/fold_change), numpy.log2(fold_change))
    for cl1, cl2 in itertools.combinations(cell_lines, 2):
        df1 = data[cl1]['median']
        df2 = data[cl2]['median']
        antibodies = df1.columns[3:]
        for ab in antibodies:
            cell_line_var = df1.loc[(df1[ab] < liml) & (df2[ab] > limu)]
