# GRAPE
GRAPE: Genetic interaction Regression Analysis of Pairwise Effects is a tool for analyzing genetic interaction scores from dual-gene CRISPR multiplex screens using regression modeling.

### Dependencies
To run this software, `Python>=3.10.1` is required. 

Software dependencies:
- `numpy==1.26.2`
- `pandas==2.1.4`
- `scipy==1.11.4`
- `statsmodels==0.14.1`
- `scikit-learn==1.3.2`

## Installation

1. Clone the repository into a directory on your computer and `cd` into the root directory of the `GRAPE` package.

```zsh
git clone https://github.com/hart-lab/GRAPE.git
cd GRAPE
```

2. Optional but recommended, create and start a virtual environment for the tool with:

```zsh
python3 -m venv grenv
source grenv/bin/activate
```


3. Install using the `setup.py` file with:

```zsh
pip3 install .
```

4. Verify installation:

```zsh
grape -h
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You should see the following:
```zsh
usage: grape [-h] [--version] -i INPUT_FILEPATH -o OUTPUT_DIRECTORY -c CONTROL_COLUMNS [CONTROL_COLUMNS ...] -t TARGET_GENE_FILE [-p OUTPUT_PREFIX] [--min-reads MIN_READS] [--pseudocount PSEUDOCOUNT] [--target-columns TARGET_COLUMNS [TARGET_COLUMNS ...]][--no-mean-replicates] [--no-groupby-targets] [--nonessential-gene-file NONESSENTIAL_GENE_FILE] [--query-gene-file QUERY_GENE_FILE] [--genepair-del GENEPAIR_DEL] [--fit-intercept] [--half-window-size HALF_WINDOW_SIZE] [--monotone-filter]

Run GRAPE

options:
  -h, --help            show this help message and exit
  --version             print version and exit
  ```

## Usage Example
```zsh
grape \
  -i path/to/read_count_file.csv \
  -o output_directory/ \
  -c T0_R1 T0_R2 \
  -t path/to/target_gene_list.txt \
  --target-columns T18_R1 T18_R2 \
  -p T18
```
## Input File Formats
### 1. Read Count file
The read count file is a **tab-delimited text file** containing raw counts for each gRNA across all replicates.
- **gRNA**: unique gRNA identifier (often `GENE1_guideindex, GENE2_guideindex, GENE1_GENE2_guideindex`).  
- **GENE**: target genes and gene pairs for the gRNA. Dual-gene constructs are joined using the delimiter specified by `--genepair-del` (default is `_`).  
- **Sample columns**: read counts for control and experimental replicates (column names must match arguments passed with `-c` and `--target-columns`).  

**Example:**
| gRNA          | GENE         | T0_R1 |T0_R2 |T18_R1|T18_R2|
| ------------- | ------------ | -----:| ----:| ----:| ----:|
| MAPK1_1       | MAPK1        | 251 | 364 | 20 | 19 |
| MAPK3_1       | MAPK3        | 445 | 724 | 85 | 31 |
| MAPK1_MAPK3_1 | MAPK1_MAPK3  | 218 | 112 | 27 | 11 |
...

### 2. Target Gene file
A plain text file with **one gene per line**, listing the target genes to include in the regression analysis (exclude controls).

**Example:** <br><br>
&nbsp;&nbsp;&nbsp;MAPK1 <br>
&nbsp;&nbsp;&nbsp;MAPK3 <BR>
&nbsp;&nbsp;&nbsp;ATM <br>
&nbsp;&nbsp;&nbsp;BLM <br>
&nbsp;&nbsp;&nbsp;...
<br><br>

## Required Arguments
| Argument                   | Description                                                   |
| -------------------------- | ------------------------------------------------------------- |
| `-i`, `--input-filepath`   | Input read count file                                         |
| `-o`, `--output-directory` | Output directory path                                         |
| `-c`, `--control-columns`  | Space-separated list of control columns (e.g., T0 replicates) |
| `-t`, `--target-gene-file` | Path to target gene list file (excluding control genes)       |

## Optional Arguments
| Argument                   | Description                                          | Default |
| -------------------------- | ---------------------------------------------------- |:-------:|
| `-p`, `--output-prefix`    | Prefix for output files                              | `None` |
| `--min-reads`              | Minimum read count threshold                         | `0` |
| `--pseudocount`            | Pseudocount to avoid division by zero                | `1` |
| `--target-columns`         | Space-separated list of target columns to average    | `None` |
| `--no-mean-replicates`     | Disable averaging across replicates                  | `False` |
| `--no-groupby-targets`     | Disable grouping by target gene                      | `False` |
| `--nonessential-gene-file` | Path to nonessential/reference gene list for mode-centering. If not provided, mode-centering is performed using the full fold-change distribution.      | `None` |
| `--query-gene-file`        | Path to query gene list                              | `None` |
| `--genepair-del`           | Delimiter for gene pairs                             | `"_"` |
| `--fit-intercept`          | Fit intercept in regression                          |`False`|
| `--half-window-size`       | Half window size for local variance, If set to 0, a global variance is calculated instead of a local one. The half-window size must NOT exceed the total number of pairwise constructs. |`500`|
| `--monotone-filter`        | Apply monotonic filter to local std deviations       |`False`|


## Output Files 
After a successful run, GRAPE produces the following files in the specified output directory (with the user-defined prefix `-p` if provided):

### 1. `grape_pairs_<prefix>.txt`
Contains regression results and genetic interaction (GI) Z-scores for gene pairs.  
- GENE_PAIR: gene pair names.  
- fc_obs: observed mode-centered fold change for the gene pair. 
- fc_exp: expected fold change predicted by the regression model for the gene pair  (sum of single-gene coefficients for gene1 and gene2).
- GI_raw: unnormalized GI score, computed as `fc_obs – fc_exp`.
- g1_fc: single-gene fold change for the first gene in the pair.
- g2_fc: single-gene fold change for the second gene in the pair.
- dLFC: observed pairwise fold change – (observed fold change of gene1 + gene2).
- local_std: local standard deviation estimated from a sliding window.
- GI_Zscore: normalized genetic interaction score; negative values indicate synthetic genetic interactions, positive values indicate suppressing interactions.
- Pval_synth: p-value for detecting synthetic (negative) interactions.
- Padj_synth: Benjamini–Hochberg adjusted p-value for synthetic interactions.
- Pval_supp: p-value for detecting suppressor (positive) interactions.
- Padj_supp: Benjamini–Hochberg adjusted p-value for suppressor interactions.

### 2. `grape_singles_<prefix>.txt`
Contains regression outputs for **single-gene effects**.  
- index: single gene names.  
- fc_obs: observed mode-centered fold change for the gene.
- fc_exp: expected fold change predicted by the regression model for the gene. (single-gene coefficient beta)

### 3. `modecenter_meanfc_<prefix>.txt`
Contains the **mode-centered mean fold-change values** used as input for regression analysis.  
- GENE: all genes and gene pairs.  
- meanFC: gene-level fold-change, averaged across all replicates. This value is mode-centered based on the method specified by the user. By default, the mode is calculated from the full FC distribution. If a `--nonessential-gene-file` is provided, the mode is calculated from the fold-change distribution of non-essential genes.
