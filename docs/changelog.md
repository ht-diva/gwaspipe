## [0.4.0] - 2026-03-18

### 🚀 Features

- Upgrade Python environment to 3.11 for gwaslab v4.0.2 compatibility, update docs and scripts
- Upgrade Python environment to 3.12
- Fill the MLOG10P column in the sumstats data using P column, restore v3 functionality

### 🐛 Bug Fixes

- Sorting imports with Ruff
- Add the GWASLab compatible logger

### 💼 Other

- Add external module for order_alleles functionality, rearrange dependencies, advance gwaslab version
- Some fixes to bypass a gwaslab buggy function

 - move forward to gwaslab 3.4.49
 - add a custom infer_build function
 - add vchange_status function from 3.6.16,
 - add version argument
- Improve the documentation
- Make pylint happy
- Improve tests
- Minor
- Merge pull request #17 from ht-diva/oa_refactor_new

Order alleles refactor

### 🚜 Refactor

- Modularize codebase by reorganizing order_alleles and utils

### 🧪 Testing

- Update test configuration

### ⚙️ Miscellaneous Tasks

- Update Python version to 3.11 for gwaslab v4.0.2 compatibility, update Dockerfile and workflow configurations
- Update test
- Ruff, sorting imports with pre-commit
- Add git-cliff config
## [0.3.1] - 2026-03-06

### 💼 Other

- Enhance README with detailed GWASPipe description

Expanded the description of GWASPipe to highlight its features and functionalities, including command-line interface and YAML configuration.
- Add to_parquet and check_ambiguous_snps
## [0.3.0] - 2025-12-12

### 🚀 Features

- Add plink-pvar format

### 💼 Other

- Read vcf study_label
- Add regenie_gene fornat for ukb-ppp
- Merge pull request #14 from ht-diva/regenie_gene

add regenie_gene format for ukb-ppp
- Add af for vcf conversion & gwascatalog_hm format
- Add info for metal_het to vcf
- Add format header for metal_het to vcf
- Change type of metal_het direction in vcf header
- Snp mapping after bcftools liftover
- Add gwascatalog custom formats
- Add gtex format
- Update snp_mapping for different SNPID formats, close #16
- Speed-up snp_mapping, close #16
- Add GENEID col
- Add GENEID to gwaslab format
- Add literature_rev format
- Add FLIPPED in snp_mapping
- Update SNPID gwaslab fun
- Update FLIPPED column
- Bump version
## [0.2.2] - 2025-05-07

### 💼 Other

- Add snp_mapping for finngen, max input decimals, decode format
- Merge pull request #13 from ht-diva/dev_finngen_float

add snp_mapping for finngen, max input decimals, decode format
- Bump version
## [0.2.1] - 2025-04-16

### 💼 Other

- Update docker-publish action
- Update README
- Update README
- Update docker playbook
- Add finngen file format
- Merge pull request #8 from ht-diva/finngen

add finngen file format
- Bump version
## [0.2.0] - 2025-01-11

### 🐛 Bug Fixes

- Fix path adding workspace_path variable
- Fix issue #4
- Fix array of arrays in yaml format , close #4
- Fix typo in config
- Fix typo
- Fix typo in the dependencies

### 💼 Other

- First stub
- Add gitignore
- Update gitignore
- Add steps and links
- Reformat
- Add python environment stuff
- Add README
- Update README
- Return two dicts for each step and add workspace folder
- Write inflation factors report as tsv file
- Update gwaspipe.py

added manhattan and qq plots step
- Merge pull request #1 from MichCM/patch-1

Update gwaspipe.py
- Update config.yml

added parameters for qq plots and manhattan plots
- Merge pull request #2 from MichCM/main

Update config.yml
- Add option to use a custom formatbook file
- Create glist-hg19.txt

created data folder and added gene list
- Create cistrans_tagger.py
- Update cistrans_tagger.py
- Update gwaspipe.py
- Update cistrans_tagger.py

minor fixes
- Update config.yml

added parameters for cistrans annotation
- Update gwaspipe.py

fixed save path for cistrans annotation (single table)
- Update config.yml
- Update cistrans_tagger.py

fixed save path with subscript for cojo or clump files
- Merge pull request #3 from MichCM/main

cis-trans annotation step
- Add custom formatbook
- Update dependencies for cistrans-tagger
- Add missing dependency
- Add newline
-  add mscorefonts to fix Matplotlib not finding basic fonts
- Update environment
- Typo
- Generalize
- Add log message
- Update README
- Add option to create workspce subfolders and report_min_pvalues step
- Specialize the vcf writing step
- Update cistrans_tagger.py

updates for the pipeline to work with minus log 10 pvalues
- Update config.yml

updates for -log10p
- Merge pull request #5 from MichCM/main

updating to work with log10P
- Mask on array to filter out study name
- Cli option to specify output path; It takes precedence over the config value
- Update config.yml
- Define a general filename mask
- Expand and add format_col_order to gwaslab format
- Add a flag for flipping alleles task and add a step to write harmonization summary
- Bump to gwaslab 3.4.43 + cache support
- Merge pull request #6 from sup3rgiu/main

Bump to gwaslab 3.4.43 + cache support
- Update the environment
- During the harmonization step, call single functions instead of the all-in-one
- Reorder the columns to be aligned with the tabix script
- Execute the flip_allele_stats as in the all-in-one function order
- Skip INFO column in regenie format as we do not trust it
- Add alphabetical sorting
- Reimplement report_min_pvalue
- Add new formats; add snp_mapping step
- Refine the configuration for processing metal heterogeneity outputs
- Update the metal_het format
- Use extreme_values to fill data
- Add pre-commit
- Add unit tests, reformat and fix minors
- Add a flag for preserving previous IDs
- Extend the gwaslab format
- Improve snp_mapping
- Specify matplotlib version

to avoid AttributeError: module 'matplotlib.cm' has no attribute 'get_cmap'
when using latest Matplotlib 3.9.0
- Update formatbook.json

adding a decode format
- Update formatbook.json
- Update formatbook.json
- Add gwaslab as input file format option
- Make a container
- Update python version requirements
- Update Dockerfile
- Update Dockerfile
- Add vcf file format
- Use short sha for the image tag
- Add study_label option for VCF input files
- Remove cis_trans_tagger
- Bump version
- Update github action
- Remove environment_dev.yml
- Merge pull request #7 from ht-diva/docker

Poetry, Docker and other stuff

### 🚜 Refactor

- Refactor
- Refactor how the parameters from the config file are handled
- Refactor and add write_pickle step
