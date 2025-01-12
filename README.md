# GWASPipe

A Python tool, based on [GWASLab](https://cloufield.github.io/gwaslab/), for assembling a computational pipeline to standardise, QC, harmonise, convert and plot summary statistics in genetic association studies (GWAS).

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

This will install snakemake into an isolated software environment

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

*   `-c` or `--config`: Path to the configuration file
*   `-i` or `--input`: Path to the summary statistics file
*   `-f` or `--format`: Format of the input file (vcf, gwaslab, regenie, fastgwa, ldsc, fuma, pickle, metal_het)
*   `-o` or `--output`: Path where results should be saved
*   `-s` or `--input_file_separator`: Input file separator
*   `-q` or `--quiet`: Set log verbosity
*   `--study_label`: Input study label, valid only for VCF files
*   `--pid`: Preserve ID

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
