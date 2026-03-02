# Getting Started with GWASPipe

This guide provides step-by-step instructions for using GWASPipe to process and harmonize GWAS summary statistics.

## Basic Usage

### Command Line Interface

```bash
gwaspipe -c config_file.yml -i input_file.tsv.gz -f input_format -o output_directory/
```

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `-c, --config_file` | YAML configuration file | `config_sumstats_harmonization.yml` |
| `-i, --input_file` | Input summary statistics file | `input_data.tsv.gz` |
| `-f, --input_file_format` | Input file format | `regenie`, `plink_pvar`, `vcf` |
| `-o, --output` | Output directory | `results/` |

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--study_label` | Study label for VCF files | None |
| `--pid` | Preserve original IDs | False |
| `--bcfliftover` | Input from BCFtools liftover | False |
| `--quiet` | Reduce log verbosity | False |

## Example Workflow

### 1. Prepare Configuration

Use the provided example configuration:

```bash
cp examples/config_sumstats_harmonization.yml my_config.yml
```

### 2. Run Processing

```bash
gwaspipe \
  -c my_config.yml \
  -i my_study.tsv.gz \
  -f regenie \
  -o my_results/ \
  --pid
```

### 3. Examine Output

The tool creates this directory structure:

```
my_results/
├── gwaspipe.log                  # Main execution log
└── outputs/
    └── my_study/                # Processed data
        ├── my_study.gwaslab.log  # Processing log
        └── my_study.gwaslab.tsv.gz  # Harmonized data
```

## Configuration File Structure

GWASPipe configuration files use YAML format with these key sections:

### Run Sequence

```yaml
run_sequence: !!omap
  - 1: 'basic_check'
  - 2: 'infer_build'
  - 3: 'fill_data'
  - 4: 'harmonize'
  - 5: 'sort_alphabetically'
  - 6: 'write_tsv'
```

### Step Configuration

Each processing step has parameters:

```yaml
steps:
  basic_check:
    params:
      run: True
    gl_params:
      threads: 4
      normalize: True

  infer_build:
    params:
      run: True
```

### Common Parameters

```yaml
# Shared parameters
n_cores: &cores 4

# Output settings
root_path: "results"
log_filename: "gwaspipe.log"
```

## Advanced Configuration

### Allele Ordering

Configure the `order_alleles` step in your workflow:

```yaml
order_alleles:
  params:
    run: True
  gl_params:
    mode: "v"  # "v" for vectorized, "p" for parallel
    n_cores: *cores
    format_snpid: True
```

### Genome Build Inference

```yaml
infer_build:
  params:
    run: True
    build: "19"  # Default build (19 or 38)
```

### Getting Help

```bash
# Show all available options
gwaspipe --help

# Check version
gwaspipe --version
```

## Next Steps

- Explore the [example configurations](examples/) for different use cases
- Review the [GWASLab documentation](https://cloufield.github.io/gwaslab/) for underlying algorithms
