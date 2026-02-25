# Getting Started

To run GWASPipe, use the following command:

## Basic Command Example:
```bash
gwaspipe -c examples/config_sumstats_harmonization.yml -i examples/input_data.tsv.gz -f regenie -o results/
```

### Usage:
```
Usage: gwaspipe [OPTIONS]

Options:
  -c, --config_file TEXT          Configuration file path  [required]
  -i, --input_file TEXT           Input file path  [required]
  -b, --formatbook_file TEXT      Formatbook file path
  -f, --input_file_format [plink_pvar|literature_rev|gtex|gwascatalog_hm_custom|ssf_custom|finngen|vcf|decode|gwaslab|regenie|regenie_gene|fastgwa|ldsc|fuma|pickle|metal_het|auto]
                                  Input file format  [required]
  -o, --output TEXT               Path where results should be saved (default=results/)
  -s, --input_file_separator TEXT
                                  Input file separator
  --study_label TEXT              Input study label, valid only for VCF files
  -q, --quiet                     Set log verbosity
  --pid                           Preserve ID
  --bcfliftover                   Input from BCFtools liftover
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```

## Expected Output

When running GWASPipe with the example configuration, the following output structure will be created in the specified output directory:

```
results/
├── gwaspipe.log                  # Main log file containing execution details
└── outputs/
    └── input_data/               # Output directory for processed data
        ├── input_data.gwaslab.log # Log file specific to the processing steps
        └── input_data.gwaslab.tsv.gz # Compressed TSV file with harmonized data
```

The main output file (`input_data.gwaslab.tsv.gz`) contains the processed and harmonized GWAS summary statistics in a standardized format. The log files provide detailed information about the processing steps and any issues encountered.


## GWASPipe Configuration File Example

The [config_sumstats_harmonization.yml](config_sumstats_harmonization.yml) configuration file is used by GWASPipe to process and harmonize summary statistics.

Each task has the following default parameter:
* `run`: indicates that this step should be performed if set to `True`, not otherwise.

 The ordered dictionary, at the top of the file, can be used to specify the order in which tasks should be executed.

### Configuration Steps:

Here's a breakdown of its tasks:

**Basic check**: This task checks for basic issues with the input data, such as missing values or incorrect formatting. It provides feedback on whether the input data needs to be corrected before proceeding.

**Infer build**: This step infers the build version of the input data. The `build` refers to the specific version of the genetic database used in the study (e.g., GRCh37 or GRCh38). This information is essential for downstream analysis and interpretation of results.

**Fill data**: In this task, GWASPipe fills missing values in the input data based on predefined rules. For example, it might use a reference distribution to fill missing values for certain variables.

It has the following parameters:
- `to_fill`: A list of column names to fill with data. In this case, it's set to `['MLOG10P', 'Z']`.
- `overwrite`: Set to `False`, which means any existing data in the specified columns will not be overwritten.
- `extreme`: Set to `True`, which means that extreme values (e.g., infinity) will be handled accordingly.

**Harmonize**: The harmonization step is where GWASPipe standardizes the input data according to specific formats or protocols. This process ensures that all input data conforms to the required format and can be easily integrated with other datasets.

**Sort alphabetically**: This task sorts the input data in alphabetical order based on a specified column (e.g., study ID). This is useful for maintaining consistency when combining data from multiple studies or experiments.

**Write TSV**: Finally, GWASPipe writes the harmonized and sorted data to a tab-separated values (TSV) file. The `fmt` parameter allows users to specify the desired output format, and other parameters control aspects like compression, float formatting, and more.

### Configuration File:

```yaml
---
# common parameters
n_cores: &cores 4

# This is an ordered dict that contains the order in which the steps are run.
run_sequence: !!omap
  - 1: 'basic_check'
  - 2: 'infer_build'
  - 3: 'fill_data'
  - 4: 'harmonize'
  - 5: 'sort_alphabetically'
  - 6: 'write_tsv'

steps:
  basic_check: # see https://cloufield.github.io/gwaslab/Standardization/
    params:
      run: True
    gl_params:
      n_cores: *cores
      normalize: True
  infer_build: # see https://cloufield.github.io/gwaslab/InferBuild/
    params:
      run: True
  fill_data: # see https://cloufield.github.io/gwaslab/Conversion/#fill_data
    params:
      run: True
    gl_params:
      to_fill: ['MLOG10P', 'Z']
      overwrite: False
      extreme: True
  harmonize: # see https://cloufield.github.io/gwaslab/Harmonization/
    params:
      run: True
    gl_params:
      basic_check: False
      n_cores: *cores
  sort_alphabetically:
    params:
      run: True
    gl_params:
      n_cores: *cores
  write_tsv:
    params:
      run: True
      workspace: "outputs"
      workspace_subfolder: True
    gl_params:
      fmt: "gwaslab"
      xymt_number: True
      to_csvargs:
        compression:
          method: "gzip"
          compresslevel: 1
          mtime: 1

# Filename transformation, e.g.: seq.10000.28.fastGWA
# Use the filename mask to extract elements from the input filename to be used in the output filename.
# If filename_mask is not provided, the input filename without the extension will be used as the output filename.
# filename_mask: [ True, True, True, False]
# filename_sep: '.'

# IO
root_path: "results"
log_filename: "gwaspipe.log"
# formatbook_path: "path/to/custom_formatbook.json"     # fill this var if you want to provide a custom formatbook
```
