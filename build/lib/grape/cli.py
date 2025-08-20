import argparse
import sys

from grape.core.run import run
from grape.utils.defaults import *
from grape.utils.version import __version__


def __main__():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Run GRAPE",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='print version and exit', action='version',
                        version=f'%(prog)s {__version__}')
    
    required_group = parser.add_argument_group('Required Arguments')
    required_group.add_argument('-i', '--input-filepath', type=str, required = True, 
                                help="Input read count file path")
    required_group.add_argument('-o', '--output-directory', type=str, required = True, help="Output " \
                                "directory path")
    required_group.add_argument('-c', '--control-columns', required=True, help='space-delimited list' \
                                'of ints or strings of control columns (T0)', nargs='+')
    required_group.add_argument('-t', '--target-gene-file', type=str, required=True, help='Path to ' \
                                'target gene list file (do not include control genes)')
    
    optional_group = parser.add_argument_group('Optional Arguments')
    optional_group.add_argument('-p', '--output-prefix', type=str, required=False, help='Prefix for ' \
                                'output files', default=DEFAULT_OUTPUT_PREFIX)
    optional_group.add_argument('--min-reads', type=int, required=False, help='Minimum read count ' \
                                'threshold for filtering samples', default=DEFAULT_MIN_READS)
    optional_group.add_argument('--pseudocount', type=int, required=False, help='Pseudocount to avoid' \
                                'devision by zero', default=DEFAULT_PSEUDOCOUNT)
    optional_group.add_argument('--target-columns', required=False, help='Space-delimited list of ' \
                                'target column names to average', nargs='+', 
                                default=DEFAULT_TARGET_COLUMNS)
    optional_group.add_argument('--no-mean-replicates', action='store_true', help='Disable averaging' \
                                'across replicate columns', default=DEFAULT_NO_MEAN_REPLICATES)
    optional_group.add_argument('--no-groupby-targets', action='store_true', help='Disable grouping ' \
                                'by target gene', default=DEFAULT_NO_GROUPBY_TARGETS)
    optional_group.add_argument('--nonessential-gene-file', type=str, help='Path to nonessential/' \
                                'reference gene list', default=DEFAULT_NONESSENTIAL_GENE_FILE)
    optional_group.add_argument('--query-gene-file', type=str, required=False, help='Path to query ' \
                                'gene list', default=DEFAULT_QUERY_GENE_FILE)
    optional_group.add_argument('--genepair-del', type=str, required=False, help = 'Delimiter used ' \
                                'to separate gene pairs in index names.', default=DEFAULT_GENEPAIR_DEL)
    optional_group.add_argument('--fit-intercept', action='store_true', help = 'Whether to fit the ' \
                                'intercept in regression', default=DEFAULT_FIT_INTERCEPT)
    optional_group.add_argument('--half-window-size', type=int, help='Half window size for local ' \
                                'variance', default=DEFAULT_HALF_WINDOW_SIZE)
    optional_group.add_argument('--monotone-filter', action='store_true', help='Apply monotonic ' \
                                'filter to local std devs', default=DEFAULT_MONOTONE_FILTER)

    arguments = parser.parse_args()
    run(arguments)


if __name__ == '__main__':
    sys.exit(__main__())
