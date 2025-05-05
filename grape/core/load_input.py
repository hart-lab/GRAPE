"""
This file handles the reading of input files.
"""

from typing import List

import pandas as pd


def load_readcount_matrix(filepath: str, index_column: int = 0, 
                          delimiter: str = '\t') -> pd.DataFrame:
    """
    Load a read count file

    Parameters: 
        filepath: str
            The path to the file to be loaded
        index_column: int
            The column to use as an index (contains target IDs). Default 0.
        delimiter: str
            tab or comma delimited file. Default '\t'
    
    Returns:
        reads_df: pd.DataFrame
            A dataframe containing the readcount matrix. Index is unique ID, first column is target.
    """  
    reads_df = pd.read_csv(filepath, index_col=index_column, sep = delimiter)
    return reads_df

def load_genelist(filepath: str) -> List[str]:
    """
    Expects a .txt file with one gene per line

    Parameters:
        filepath: str
            The path to the genelist
 
    Returns:
        List of strings
    """
    genes = []
    with open(filepath, 'r', encoding='utf-8') as infile:
        for line in infile:
            gene = line.strip('\n')
            genes.append(gene)
	
    return genes