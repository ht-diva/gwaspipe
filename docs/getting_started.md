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
  - 5: 'canonicalize_effect_alleles'
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

### Assembly declaration and provenance

Reference-aware steps (`harmonize`, `liftover`, and VCF export) require a preceding
`infer_build` step with a declared input assembly. `infer_build` compares HapMap3
coordinate matches for GRCh37 and GRCh38, then records its decision in GWASLab
metadata and in a `.provenance.json` sidecar next to each output. Declare the
versions of every reference resource used by the workflow.

```yaml
genome_assembly: GRCh38
reference_resources:
  hapmap3_coordinates: "GWASLab 4.0.2 bundled HapMap3 tables"
  reference_fasta: "GRCh38, release 109"
  allele_frequency_panel: "1000 Genomes Phase 3, GRCh38"
assembly_validation:
  min_hapmap3_matches: 10000
  allow_override: false
```

An override is permitted only when it is explicitly justified and is recorded in
the provenance:

```yaml
assembly_validation:
  allow_override: true
  override_reason: "Legacy study has fewer than 10,000 HapMap3 markers; validated against its manifest."
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

Configure the `canonicalize_effect_alleles` step in your workflow. Its parameters are
forwarded directly to `SumstatsManager.order_alleles()`:

```yaml
canonicalize_effect_alleles:
  params:
    run: True
  gl_params:
    mode: "v"  # "v" for vectorized, "p" for parallel
    n_cores: *cores
    format_snpid: True
```

`sort_alphabetically` remains available as a deprecated compatibility alias.

### Genome Build Inference

```yaml
infer_build:
  params:
    run: True
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
