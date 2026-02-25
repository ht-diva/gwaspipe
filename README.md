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

## Documentation

You can use one of the following ways for installing Gwaspipe.

### Installation via Conda/Mamba

1.  Clone the repository using Git: `git clone https://github.com/ht-diva/gwaspipe.git`
2.  Create a conda environment: `conda env create -f environment.yml`
3.  Activate it: `conda activate gwaspipe`
4.  Install gwaspipe: `make install`

This will install gwaspipe into an isolated software environment

5.  Run the tool: `gwaspipe --help`

### Docker image

The Linux OS image is available from the github packages repository:
[Docker image](https://github.com/ht-diva/gwaspipe/pkgs/container/gwaspipe)


### Getting Started

To get started with GWASPipe have a look at the [getting started guide](docs/getting_started.md)

### Contributing

We welcome contributions to GWASPipe! If you'd like to contribute, please follow these guidelines:

1.  Fork the repository on GitHub: `https://github.com/your-username/gwaspipe.git`
2.  Create a new branch for your changes: `git checkout -b feature/your-feature-name`
3.  Make your changes and commit them: `git add . && git commit -m "Your commit message"`
4.  Push your changes to the remote repository: `git push origin feature/your-feature-name`
5.  Open a pull request on GitHub to propose your changes

We appreciate your contributions and look forward to seeing GWASPipe grow!
