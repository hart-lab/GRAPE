"""
This file calculates fold change from read count file
"""


from typing import List, Union, Optional

import numpy as np
import pandas as pd
import scipy.stats as stats


def get_foldchange_matrix(reads_df: pd.DataFrame, control_columns: Union[List[str], List[int]], 
						  min_reads: int = 0, pseudocount: int = 1) -> pd.DataFrame:
	"""
    Given a dataframe of raw read counts,
    1. filter for T0 min read counts.
    2. calculate log2 ratio of each sample relative to the mean of control columns.
    
    Parameters:
	    reads_df : pd.DataFrame
            DataFrame containing raw read counts, with index as unique IDs, columns as samples,
			and first column is the target-gene column
        control_columns : List[str] or List[int]
            Either a list of control sample labels (column names) or a list of their column indices
        min_reads : int, optional
            Minimum read count threshold for filtering samples #TODO(Juihsuan) add later. Default is 0.
        pseudocount : int, optional
            Small value added to read counts to avoid division by zero and log transformations of 
			zero. Default is 1.

    Returns:
        pd.DataFrame
            DataFrame with log2 fold-change values, where index is the same as `reads_df` 
            and columns correspond to target samples
    """
    # control_columns to actual column names
	try:
		column_list = list(map(int, control_columns))
		control_column_labels = reads_df.columns.values[column_list]
	except ValueError:
		control_column_labels = control_columns

	print("[INFO] Using controls: " + ",".join(map(str, control_column_labels)))

	target_column_labels = [x for x in reads_df.columns.values if x not in control_column_labels]
	# target_column_labels[0] should still be the target gene

	ctrl_sum = reads_df[control_column_labels].sum(axis=1)
	fc_df = pd.DataFrame(index=reads_df.index.values, columns=target_column_labels, data=0.)
	fc_df[target_column_labels[0]] = reads_df[target_column_labels[0]] 
	for col_name in target_column_labels[1:]:
		fc_df[col_name] = np.log2(((reads_df[col_name].values + pseudocount)/
							        sum(reads_df[col_name].values)) /
		                           ((ctrl_sum + pseudocount) / sum(ctrl_sum)))

	return fc_df

def get_mean_foldchange(fc_df: pd.DataFrame, target_columns: Optional[List[str]] = None, 
						no_mean_replicates: bool = False, no_groupby_targets: bool = False) -> pd.DataFrame:
	"""
    Given a dataframe of fold changes, return the mean across replicates and/or the mean
    across guides targeting the same gene(s)
    
    Parameters:
	    fc_df : pd.DataFrame
            DataFrame of log2 fold-change values, where the first column is the target gene(s).
			Subsequent columns are the replicates to be averaged.
        target_columns : List[str], optional
            List of column labels indicating replicates to be averaged. If None, all columns except 
			the first are used.
        no_mean_replicates : bool
            Whether to average across replicate columns. Default is False.
        no_groupby_targets : bool, optional
            Whether to compute the mean fold-change by target gene(s). Default is False.

    Returns:
        pd.DataFrame
            DataFrame with mean fold-change values, optionally grouped by target genes
    """
	if target_columns:
		# use target columns
		target_columns = target_columns
	else:
		# use all columns
		target_columns = fc_df.columns.values

	if no_mean_replicates:
		outcols = fc_df.columns.values
	else:
		outcols = [fc_df.columns.values[0], 'meanFC']

	mean_fc_df = pd.DataFrame( index = fc_df.index.values, columns=outcols, data=0.)
	mean_fc_df[outcols[0]] = fc_df[outcols[0]]

	if no_mean_replicates:
		mean_fc_df[outcols[1:]] = fc_df[fc_df.columns.values[1:]]
	else:
		mean_fc_df[outcols[1]] = fc_df[target_columns].mean(1)

	if not no_groupby_targets:
		mean_fc_df = mean_fc_df.groupby(outcols[0]).mean()

	return mean_fc_df

def mode_center(mean_fc_df: pd.DataFrame) -> pd.DataFrame:
	"""
	Assumes a polished fold change df where the index is the target gene(s). Normalize fold-change 
	mode to zero across entire distribution.
	
	Parameters:
        mean_fc_df: pd.DataFrame
            DataFrame containing fold-change values, where the index represents target genes
			
	Returns:
        pd.DataFrame
            Mode-centered fold-change values
    """
	xx = np.linspace(-5, 4, 901)
	fc_col = mean_fc_df.columns[-1]
	kx = stats.gaussian_kde(mean_fc_df[fc_col].astype(float))
	mode_x = xx[np.argmax(kx.evaluate(xx))]

	mean_fc_df[fc_col] = mean_fc_df[fc_col] - mode_x

	return mean_fc_df

def mode_center_vs_reference_genes(mean_fc_df: pd.DataFrame, noness_genes: List[str]) \
	-> pd.DataFrame:
	"""
	Assumes a polished fold change df where the index is the target gene(s). Normalize fold-change 
	values using the median of a reference set of non-essential genes.
	
	Parameters:
        mean_fc_df: pd.DataFrame
            DataFrame containing fold-change values, where the index represents target genes
		noness_genes: List[str]
            List of non-essential genes
			
	Reeturns:
        pd.DataFrame
            Fold-change values normalized relative to the median of reference genes
    """
	nonidx = [x for x in mean_fc_df.index.values if x in noness_genes]
	modecenter_fc_df = mean_fc_df - mean_fc_df.loc[nonidx].median()

	return modecenter_fc_df
