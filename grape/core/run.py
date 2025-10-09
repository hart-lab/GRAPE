"""
This file holds the workflow for the GRAPE package
"""


from argparse import Namespace

from grape.core.load_input import load_genelist, load_readcount_matrix
from grape.core.foldchange_generator import *
from grape.core.regression import *
from grape.core.zscore_generator import *


def run(args: Namespace) -> None:

    # 01: load input read count file
    reads = load_readcount_matrix(args.input_filepath)

    # 02: generate fold change input
    raw_fc = get_foldchange_matrix(reads, args.control_columns, args.min_reads, args.pseudocount)
    fc = get_mean_foldchange(raw_fc, args.target_columns, args.no_mean_replicates, 
                             args.no_groupby_targets)
    if args.nonessential_gene_file is None:
        modecenter_meanfc = mode_center(fc)
    else:
        noness = load_genelist(args.nonessential_gene_file)
        modecenter_meanfc = mode_center_vs_reference_genes(fc, noness)
    
    # 03: load target gene list
    target_gene_list = load_genelist(args.target_gene_file)
    if args.query_gene_file is not None:
        query_genes = load_genelist(args.query_gene_file)
        target_gene_list = list(set(target_gene_list + query_genes))

    # 04: run regression
    predictor, obs = make_predictor_matrix(modecenter_meanfc, target_gene_list, args.genepair_del)
    pairs, singles, model = do_regression(predictor, obs, args.fit_intercept, args.genepair_del)

    # 05: apply dynamic range filter and re-run regression
    remove_me = dynamic_range_filter(pairs)
    print(f'[INFO] Dynamic range filter removed {len(remove_me)} gene pairs.')

    fc_filt = fc.drop(remove_me.index.values, axis=0, errors='ignore')

    if args.nonessential_gene_file is None:
        modecenter_meanfc_filt = mode_center(fc_filt)
    else:
        noness = load_genelist(args.nonessential_gene_file)
        modecenter_meanfc_filt = mode_center_vs_reference_genes(fc_filt, noness)

    predictor, obs = make_predictor_matrix(modecenter_meanfc_filt, target_gene_list, args.genepair_del)
    pairs, singles, model = do_regression(predictor, obs, args.fit_intercept, args.genepair_del)
    print(
        f"[INFO] Regression R²: {model['Rsq']:.3f}, Intercept: {model['Intercept']}")
    
    # 06: calculate GI Zscore
    pairs_localZ = get_zscore(pairs, args.half_window_size, args.monotone_filter)
    
    # 07: save outputs
    output = pairs_localZ.copy()
    float_format = '{:.3f}'
    sci_format = '{:.2e}'
    output.index.name='GENE_PAIR'
    
    for col in output.columns[:8]:
        output[col] = output[col].map(lambda x: float_format.format(x))

    for col in output.columns[8:]:
        output[col] = output[col].map(lambda x: sci_format.format(x))

    output_path = args.output_directory.rstrip('/') + '/'
    prefix = f"_{args.output_prefix.rstrip('_')}" if args.output_prefix else ""
    
    output.to_csv(output_path + f'grape_pairs{prefix}.txt', sep='\t')
    
    singles_fmt = singles.copy()
    singles_fmt = singles_fmt.map(lambda x: float_format.format(x))
    singles_fmt.to_csv(output_path + f'grape_singles{prefix}.txt', sep='\t')
    
    modecenter_meanfc_fmt = modecenter_meanfc.copy()
    modecenter_meanfc_fmt = modecenter_meanfc_fmt.map(lambda x: float_format.format(x))
    modecenter_meanfc_fmt.to_csv(output_path + f"modecenter_meanfc{prefix}.txt", sep="\t")

    print(f"[INFO] GRAPE analysis complete. Results saved to {output_path}.")

