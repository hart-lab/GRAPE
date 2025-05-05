"""
This file makes predictor matrix, run regression, and filter regression
"""


from typing import Tuple, Dict, List

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


def make_predictor_matrix(fc_df: pd.DataFrame, target_gene_list: List[str], 
						  gene_delim: str = '_') -> Tuple[pd.DataFrame, pd.DataFrame]:
	"""
	Constructs a binary predictor matrix for regression analysis.
    
    Parameters:
        fc_df : pd.DataFrame
            DataFrame containing fold-change values with index as genes and gene pairs
        target_gene_list : List[str]
            List of target genes to be included in the regression model
        gene_delim : str
            Delimiter used to separate gene pairs in index names. Default is "_"

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]
            - Predictor matrix (binary matrix where columns are target genes and rows are gene 
			  knockouts)
            - Response vector (observed fold-change values)
	"""
	predictor_matrix = pd.DataFrame(index=fc_df.index.values, columns=target_gene_list, data=0.)
	for target_array in fc_df.index.values:
		if gene_delim in target_array:
			
            # this target is a gene pair
			target_genes_in_array = target_array.split(gene_delim)
			if (len(np.intersect1d(target_genes_in_array, target_gene_list)) == \
	            len(target_genes_in_array)):
				# both/all genes in the target array are in the genetic interaction gene list. 
                # include this genepair.
				predictor_matrix.loc[target_array, target_genes_in_array] = 1
		else:
			# absence of delimiter implies single gene knockout. is the single gene in our 
            # target_gene_list?
			if (target_array in target_gene_list):
				predictor_matrix.loc[target_array, target_array] = 1
	
	# Drop rows not targeting a gene in the target list
	dropme = np.where( predictor_matrix.sum(1)==0 )[0]  
	predictor_matrix.drop( predictor_matrix.index.values[dropme], axis=0, inplace=True )

	# Drop columns (genes) with no arrays targeting that gene
	dropme = np.where( predictor_matrix.sum(0)==0 )[0]  
	predictor_matrix.drop( predictor_matrix.columns.values[dropme], axis=1, inplace=True )

	obs_vector = fc_df.loc[ predictor_matrix.index.values ]
	print('[INFO] regression matrix rows: {:5d}, cols: {:3d}'.format(predictor_matrix.shape[0], 
															predictor_matrix.shape[1]))
	return predictor_matrix, obs_vector

def do_regression(predictor_matrix: pd.DataFrame, obs_vector: pd.DataFrame, 
				  fit_intercept: bool = False, delimiter: str = '_') -> \
                  Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]:
	"""
    Calculate the regression and provide the initial GI score.

    Parameters:
        predictor_matrix : pd.DataFrame
            Binary predictor matrix generated from make_predictor_matrix
        obs_vector : pd.DataFrame
            Response vector generated from make_predictor_matrix
        fit_intercept : bool
            Whether to fit an intercept in the regression model. Default is False.
        delimiter : str
            Delimiter used to separate gene pairs in index names. Default is "_".

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]
            - DataFrame containing predicted fold-changes, genetic interaction scores, and LFC
            - DataFrame containing single-gene LFC
            - Dictionary with regression metadata (R-squared value, intercept, and model parameters)	
    """
	model = LinearRegression(fit_intercept=fit_intercept).fit(predictor_matrix.values, obs_vector)
	pred_fc = model.predict(predictor_matrix.values)
	
	pairs = pd.DataFrame(index=predictor_matrix.index.values, 
					     columns=['fc_obs','fc_exp','GI_raw','g1_fc','g2_fc','dLFC'], data=0.)
	pairs['fc_obs']  = obs_vector.values
	pairs['fc_exp'] = pred_fc.flatten()  # already an array
	
    # remove singles: columns of predictor matrix
	single_genes = predictor_matrix.columns.values
	singles = pairs.loc[single_genes,['fc_obs','fc_exp']]
	pairs.drop( single_genes, axis=0, inplace=True )

	pairs['GI_raw'] = pairs.fc_obs - pairs.fc_exp

	# now for each pair, get g1, g2, dLFC based on observed FC in singles
	for genepair in pairs.index.values:
		g1, g2 = genepair.split(delimiter)
		g1_fc = singles.loc[g1,'fc_obs']
		g2_fc = singles.loc[g2,'fc_obs']
		dLFC = pairs.loc[genepair, 'fc_obs'] - (g1_fc + g2_fc)
		pairs.loc[genepair, ['g1_fc','g2_fc','dLFC']] = g1_fc, g2_fc, dLFC

	# other data:
	metadata = {}
	metadata['Rsq'] = model.score(predictor_matrix.values, obs_vector.values)
	metadata['Intercept'] = model.intercept_
	metadata['Params']  = model.get_params()
	
	return pairs, singles, metadata

def dynamic_range_filter(regression_pairs: pd.DataFrame) -> pd.DataFrame:
	"""
	Identify target list where expected phenotype (fc_exp) is beyond the dynamic range of the assay 
	(min fc_obs). These break the regression model.

	Parameters:
	    regression_pairs: pd.DataFrame
		    'pairs' output from do_regression()
	
	Returns:
	    subset of regression_pairs dataframe whose indices should be deleted from the foldchange df,
		and run grape do_regression again.
	"""
	fc_limit = regression_pairs['fc_obs'].min()
	pairs_to_remove = regression_pairs[regression_pairs['fc_exp'] < fc_limit]
	
	return pairs_to_remove
