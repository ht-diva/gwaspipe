# GWASPipe

**GWASPipe** is a Python-based tool that streamlines the management of **genetic association study data** by automating complex workflows. The tool serves as a **computational pipeline** that enables researchers to perform essential tasks such as **quality control, standardization, and visualization** of summary statistics. Users interact with the program via a **command-line interface**, using customizable **YAML configuration files** to tailor data processing to their specific needs.
GWASPipe composes reusable processing steps (normalization, allele harmonization, filtering, QC metrics, plotting) using established libraries such as [GWASLab](https://cloufield.github.io/gwaslab/) and pandas, producing cleaned, analysis-ready summary files along with QC reports and figures—making multi-study harmonization reproducible, efficient, and less error-prone.
### Table of Contents

*   [Requirements](#requirements)
*   [Getting Started](#getting-started)
*   [Usage](#usage)
*   [Example Configuration Files](#example-configuration-files)
*   [Contributing](#contributing)

## Requirements

To use GWASPipe, you'll need:

*   Python 3.10 or higher
*   The following dependencies installed:
    *   `click` for command-line interface
    *   `cloup` for file handling and metadata management
    *   `loguru` for logging
    *   `ruamel-yaml` for YAML parsing and generation
    *   `pandas` for data manipulation and analysis
    *   `pyarrow` for high-performance data processing
    *   `numpy` for numerical computations
    *   `matplotlib` for plotting
    *   `gwaslab` for handling sumstats

see [pyproject.toml](environment.yml)

## Getting Started

You can use one of the following ways for installing Gwaspipe.

### Installation via Conda/Mamba

1.  Clone the repository using Git: `git clone https://github.com/your-username/gwaspipe.git`
2.  Create a conda environment: `conda env create -n gwaspipe -f environment.yml`
3.  Activate it: `conda activate gwaspipe`
4.  Install gwaspipe: `make install`

This will install gwaspipe into an isolated software environment

5.  Run the tool: `gwaspipe --help`

### Docker image

Gwaspipe is available also as a [Dockerfile](Dockerfile)

1.  Build the image locally: `docker build -t gwaspipe:latest .`
2.  Run it: `docker run -t -i gwaspipe:latest gwaspipe --help`

The Linux OS image is available from the github packages repository:
[Docker image](https://github.com/ht-diva/gwaspipe/pkgs/container/gwaspipe)
## Usage

GWASPipe provides a command-line interface (CLI) for easy usage. You can customize the behavior of the tool by providing configuration files in YAML format.

The CLI takes the following arguments:

* `-c` or `--config_file`: Path to the configuration file
* `-i` or `--input_file`: Path to the summary statistics file
* `-b` or `--formatbook_file`: Formatbook file path
* `-f` or `--input_file_format`: Format of the input file (plink_pvar, literature_rev, gtex, gwascatalog_hm_custom, ssf_custom, finngen, vcf, decode, gwaslab, regenie, regenie_gene, fastgwa, ldsc, fuma, pickle, metal_het, auto)
* `-o` or `--output`: Path where results should be saved
* `-s` or `--input_file_separator`: Input file separator
* `-q` or `--quiet`: Set log verbosity
* `--study_label`: Input study label, valid only for VCF files
* `--pid`: Preserve ID - When enabled, creates additional columns in the output to preserve original SNP identifiers
* `--bcfliftover`: Input from BCFtools liftover

Example configuration files are provided in the `examples` directory.

### Example Configuration Files

The following example configuration files demonstrate how to customize GWASPipe:

*   [config_harmonize_VCF_sumstats.yml](examples/config_harmonize_VCF_sumstats.yml.md): a sample configuration file for harmonizing summary statistics in [gwas-vcf](https://github.com/MRCIEU/gwas-vcf-specification) format


### Contributing

We welcome contributions to GWASPipe! If you'd like to contribute, please follow these guidelines:

1.  Fork the repository on GitHub: `https://github.com/your-username/gwaspipe.git`
2.  Create a new branch for your changes: `git checkout -b feature/your-feature-name`
3.  Make your changes and commit them: `git add . && git commit -m "Your commit message"`
4.  Push your changes to the remote repository: `git push origin feature/your-feature-name`
5.  Open a pull request on GitHub to propose your changes

We appreciate your contributions and look forward to seeing GWASPipe grow!
