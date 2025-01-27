import os
import pandas as pd

def parse_wormbase_mitab(file_path: str) -> pd.DataFrame:
    """Read a miTab TSV file of C. elegans gene interactions
     and convert it into a Pandas DataFrame.

    Parameters
    ----------
    file_path : str
        Path to TSV file that is to be read.

    Returns
    -------
    dataframe : pd.DataFrame
        DataFrame of interaction data.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    try:
        # Step 1: Locate the header line (last line starting with '#')
        with open(file_path, 'r') as file:
            lines = file.readlines()

        header_line = None
        header_index = None

        for i, line in enumerate(lines):
            if line.startswith('#') and not line.strip().startswith('######'):
                header_line = line.strip('#').strip()
                header_index = i

        if header_line is None:
            raise ValueError("Header line not found in the TSV file.")

        columns = header_line.split('\t')
        with open(file_path, "r") as file:
            for i, line in enumerate(file):
                if len(line.split("\t")) != len(columns):
                    print(f"Issue at line {i}: {line}")

        dataframe = pd.read_csv(file_path, sep='\t', skiprows=header_index + 1,
                                names=columns, na_values="-")
        if dataframe.columns.isnull().any():
            raise ValueError("Invalid column headers in the TSV file.")
        return dataframe
    except Exception as e:
        raise ValueError(f"Error reading TSV file: {e}")