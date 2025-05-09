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
git clone https://github.com/yourusername/GRAPE.git
cd GRAPE
```

2. Optional but recommended, create and start a virtual environment for the tool with:

```zsh
python3 -m venv grenv
source grenv/bin/activate
```

To exit the environment later, run `deactivate` in the command line. 

3. Install using the `setup.py` file with:

```zsh
pip3 install .
```

4. Verify installation:

```zsh
grape -h
```
You should see the following:
```zsh
usage: grape [-h] [--version] -i INPUT_FILEPATH -o OUTPUT_DIRECTORY -c CONTROL_COLUMNS [CONTROL_COLUMNS ...] -t TARGET_GENE_FILE
             [--min-reads MIN_READS] [--pseudocount PSEUDOCOUNT] [--target-columns TARGET_COLUMNS [TARGET_COLUMNS ...]]
             [--no-mean-replicates] [--no-groupby-targets] [--nonessential-gene-file NONESSENTIAL_GENE_FILE]
             [--query-gene-file QUERY_GENE_FILE] [--genepair-del GENEPAIR_DEL] [--fit-intercept] [--half-window-size HALF_WINDOW_SIZE]
             [--monotone-filter]

Run GRAPE

options:
  -h, --help            show this help message and exit
  --version             print version and exit
  ```

## Quick Usage Example
```zsh
grape \
  -i path/to/read_count_file.csv \
  -o output_directory/ \
  -c T0_R1 T0_R2 T0_R3 \
  -t path/to/target_gene_list.txt \
  --target-columns T18_R1 T18_R2 T18_R3
```
### Required Arguments
| Argument                   | Description                                                   |
| -------------------------- | ------------------------------------------------------------- |
| `-i`, `--input-filepath`   | Input read count file                                         |
| `-o`, `--output-directory` | Output directory path                                         |
| `-c`, `--control-columns`  | Space-separated list of control columns (e.g., T0 replicates) |
| `-t`, `--target-gene-file` | Path to target gene list file (excluding control genes)       |

### Optional Arguments
| Argument                   | Description                                          | Default |
| -------------------------- | ---------------------------------------------------- |:-------:|
| `--min-reads`              | Minimum read count threshold                         | `0` |
| `--pseudocount`            | Pseudocount to avoid division by zero                | `1` |
| `--target-columns`         | Space-separated list of target columns to average    | `None` |
| `--no-mean-replicates`     | Disable averaging across replicates                  | `False` |
| `--no-groupby-targets`     | Disable grouping by target gene                      | `False` |
| `--nonessential-gene-file` | Path to nonessential/reference gene list for mode-centering. If not provided, mode-centering is performed using the full fold-change distribution. | `None` |
| `--query-gene-file`        | Path to query gene list                              | `None` |
| `--genepair-del`           | Delimiter for gene pairs                             | `"_"` |
| `--fit-intercept`          | Fit intercept in regression       |`False`|
| `--half-window-size`       | Half window size for local variance  |`500`|
| `--monotone-filter`        | Apply monotonic filter to local std deviations       |`False`|
