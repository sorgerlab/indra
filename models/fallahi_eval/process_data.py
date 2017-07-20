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


def find_extremes(data, fold_change, save_file=None):
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
                all_extremes.append([drug, time, conc, cell_line, ab, val])
    # Sort values
    all_extremes = sorted(all_extremes, key=lambda x: abs(x[5]), reverse=True)
    # Optionally save into a CSV file
    if save_file:
        with open(save_file, 'w') as fh:
            fh.write('Drug,Time (hr),Concentration (uM),' + \
                     'CellLine,Antibody,Value\n')
            for vals in all_extremes:
                fh.write(','.join([str(v) for v in vals]) + '\n')
    return all_extremes


def find_cell_line_vars(data, fold_change, save_file=None):
    """Return conditions in which cell lines are qualitatively different."""
    liml, limu = (numpy.log2(1.0/fold_change), numpy.log2(fold_change))
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
    all_vals = sorted(all_vals, key=lambda x: abs(x[6]-x[7]), reverse=True)
    # Optionally save into a CSV file
    if save_file:
        with open(save_file, 'w') as fh:
            fh.write('Drug,Time (hr),Concentration (uM),Antibody,' + \
                     'CellLine1,CellLine2,Value1,Value2\n')
            for vals in all_vals:
                fh.write(','.join([str(v) for v in vals]) + '\n')
    return all_vals
