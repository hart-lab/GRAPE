"""
This file contains definitions to run the tool
"""

from typing import Optional, List

DEFAULT_OUTPUT_PREFIX: str = None
"""Prefix for output files"""

DEFAULT_MIN_READS: int = 0
"""Minimum read count threshold for filtering samples"""

DEFAULT_PSEUDOCOUNT: int = 1
"""Pseudocount added to avoid division by zero when computing fold changes"""

DEFAULT_TARGET_COLUMNS: Optional[List[str]] = None
"""Space-delimited list of target column names to average"""

DEFAULT_NO_MEAN_REPLICATES: bool = False
"""If True, disables averaging across replicate columns"""

DEFAULT_NO_GROUPBY_TARGETS: bool = False
"""If True, disables grouping fold changes by target gene"""

DEFAULT_NONESSENTIAL_GENE_FILE: str = None
"""Path to nonessential/reference gene list"""

DEFAULT_QUERY_GENE_FILE: str = None
"""Path to query gene list"""

DEFAULT_GENEPAIR_DEL: str = '_'
"""Default delimiter for separating gene pairs in index names"""

DEFAULT_FIT_INTERCEPT: bool = False
"""If True, fits intercept in regression models"""

DEFAULT_HALF_WINDOW_SIZE: int = 500
"""Half-window size used when computing local variance for GI Z-scores"""

DEFAULT_MONOTONE_FILTER: bool = False
"""If True, applies monotonic filtering to GI Z-score computation"""
