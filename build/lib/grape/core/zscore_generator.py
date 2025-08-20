"""
This file calculates Z score of genetic interaction and p-values, FDR for suppressing and synthetic 
interactions.
"""


import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection


def get_zscore(regression_df: pd.DataFrame, half_window_size: int = 500, 
			   monotone_filter: bool = False):
	"""
	Calculate zscore of genetic interactions through drugz variance window method

	Parameters: 
        regression_df : pd.DataFrame
            DataFrame containing genetic interaction data, output from `do_regression()`
        half_window_size : int, optional, default=500
            Half-window size for calculating local variance. If set to 0, calculates global variance
        monotone_filter : bool, optional, default=False
            If True, ensures a monotonically increasing variance.
	
	Returns:
        pd.DataFrame
            Updated DataFrame with added columns:
            - `local_std`: Local standard deviation of genetic interactions.
            - `GI_Zscore`: Z-score of genetic interactions.
            - `Pval_synth`: P-value for synthetic interactions.
            - `Padj_synth`: Adjusted p-value (FDR) for synthetic interactions.
            - `Pval_supp`: P-value for suppressing interactions.
            - `Padj_supp`: Adjusted p-value (FDR) for suppressing interactions.
	"""
	# sort the regression DF by expected fold change:
	zscore_df = regression_df.sort_values('fc_exp', ascending=False).copy()

	# Intialize output columns
	zscore_df[['local_std','GI_Zscore']] = 0.

	if (half_window_size==0):
		# if half_window_size = 0, use global instead of local Z score
		zscore_df['local_std'] = zscore_df['GI_raw'].std()
		zscore_df['GI_Zscore'] = stats.zscore( zscore_df.GI_raw.values)
	else:
		# otherwise step through the data and calculate local std
		stepsize = int( np.ceil(half_window_size / 5))
		for idx in range(half_window_size, len(zscore_df) - half_window_size, stepsize):
			
			# select the window slice
			bin_ = zscore_df.iloc[idx - half_window_size : idx + half_window_size ].copy()

			# remove outliers and calculate std
			Q1,Q3 = bin_['GI_raw'].quantile([0.25,0.75])
			IQR = Q3 - Q1
			lower_bound = Q1 - 1.5*IQR
			upper_bound = Q3 + 1.5*IQR
			bin_ = bin_[(bin_['GI_raw'] >= lower_bound) & (bin_['GI_raw'] <= upper_bound)]
			local_std = bin_['GI_raw'].std()

			# assigne value in df, 
			if (monotone_filter):
				prev_std = zscore_df.iloc[idx-1, zscore_df.columns.get_loc('local_std')]
				zscore_df.iloc[idx:idx+stepsize, zscore_df.columns.get_loc('local_std')]= \
				                                                max(local_std, prev_std)
			else:
				zscore_df.iloc[idx:idx+stepsize, zscore_df.columns.get_loc('local_std')]=local_std

		# now for the first (0..half_window_size) and last segments of the df
		zscore_df.iloc[:half_window_size, zscore_df.columns.get_loc('local_std')] = \
		                zscore_df.iloc[half_window_size, zscore_df.columns.get_loc('local_std')]
		zscore_df.iloc[-half_window_size:, zscore_df.columns.get_loc('local_std')]= \
		                zscore_df.iloc[-half_window_size-1, zscore_df.columns.get_loc('local_std')]
		
		# and calculate z score
		zscore_df['GI_Zscore'] = zscore_df['GI_raw'] / zscore_df['local_std']  # mean had better be zero
		
	# calculate synthetic and suppressing interaction p-value and FDR
	zscore_df.sort_values('GI_Zscore', ascending=True, inplace=True)
	zscore_df['Pval_synth'] = stats.norm.cdf(zscore_df.GI_Zscore)
	zscore_df['Padj_synth'] = fdrcorrection(zscore_df.Pval_synth)[1] 
    
	zscore_df.sort_values('GI_Zscore', ascending=False, inplace=True)
	zscore_df['Pval_supp'] = stats.norm.sf(zscore_df.GI_Zscore)
	zscore_df['Padj_supp'] = fdrcorrection(zscore_df.Pval_supp)[1]
    
	return zscore_df.sort_values('GI_Zscore', ascending=True)