# GWASPipe

**GWASPipe** is a Python-based computational pipeline that streamlines the management and processing of **genome-wide association study (GWAS) summary statistics**. It automates complex workflows for quality control, standardization, and visualization, making multi-study harmonization more reproducible, efficient, and less error-prone.

### Key Features

- **Modular Architecture**: Organized into reusable components (`order_alleles`, `utils`)
- **Automated Workflows**: Handles normalization, allele harmonization, filtering, and QC metrics
- **Flexible Configuration**: Uses YAML configuration files for customizable processing
- **Comprehensive Reporting**: Generates QC reports and visualizations
- **High Performance**: Leverages parallel processing and optimized algorithms

### Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Module Documentation](#module-documentation)
- [Contributing](#contributing)
- [License](#license)

## Requirements

- **Python**: 3.11 or higher
- **Dependencies**: See `pyproject.toml` for complete list
- **Key Packages**:
  - `gwaslab` - Core GWAS processing library
  - `pandas`, `numpy` - Data manipulation
  - `click`, `cloup` - Command-line interface
  - `loguru` - Advanced logging
  - `ruamel.yaml` - YAML configuration

## Installation

### Via Conda/Mamba (Recommended)

```bash
# Clone the repository
git clone https://github.com/ht-diva/gwaspipe.git
cd gwaspipe

# Create and activate environment
conda env create -f environment_docker.yml
conda activate gwaspipe

# Install package
make install
```

### Via Docker

```bash
# Pull the Docker image
docker pull ghcr.io/ht-diva/gwaspipe:latest

# Run container
docker run -v $(pwd):/data ghcr.io/ht-diva/gwaspipe gwaspipe --help
```

## Quick Start

Process GWAS summary statistics with a single command:

```bash
gwaspipe \
  -c examples/config_sumstats_harmonization.yml \
  -i examples/input_data.tsv.gz \
  -f regenie \
  -o results/
```

## Module Documentation

### Order Alleles Module

The `gwaspipe.order_alleles` module provides comprehensive allele ordering functionality:

```python
from gwaspipe.order_alleles import order_alleles
import pandas as pd
from gwaslab.g_Log import Log

# Basic usage
df = pd.DataFrame({
    'CHR': [1, 2, 3],
    'POS': [1000, 2000, 3000],
    'EA': ['A', 'T', 'C'],
    'NEA': ['T', 'A', 'G'],
    'STATUS': [9999999, 9999999, 9999999]
})

log = Log()
result = order_alleles(df, log=log)
```


## Configuration

GWASPipe uses YAML configuration files to define processing pipelines. See [Getting Started Guide](docs/getting_started.md) for detailed configuration examples.

## Contributing

We welcome contributions! Please follow these steps:

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/your-feature`
3. **Make changes** and add tests
4. **Commit changes**: `git commit -m "Add feature description"`
5. **Push branch**: `git push origin feature/your-feature`
6. **Open a Pull Request**

## License

GWASPipe is released under the [MIT License](LICENSE).
