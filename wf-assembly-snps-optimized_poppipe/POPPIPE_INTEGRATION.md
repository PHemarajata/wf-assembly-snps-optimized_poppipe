# PopPIPE-bp Integration

This document describes the integration of PopPIPE-bp outputs with the wf-assembly-snps pipeline for cluster-based SNP analysis and tree grafting.

## Overview

The PopPIPE-bp integration allows you to:

1. **Parse PopPIPE-bp outputs** - Automatically extract cluster assignments and strain information
2. **Run SNP analysis per cluster** - Perform detailed SNP analysis on each cluster separately
3. **Tree grafting** - Combine cluster-specific phylogenies using the PopPIPE tree grafting algorithm
4. **Generate comprehensive results** - Produce integrated results with accurate SNP representation across all samples

## Usage

### Basic Command

```bash
# Using TSV file with sample-to-file mapping (recommended)
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe/output \
  --input /path/to/poppipe/rfiles.txt \
  --outdir results \
  --snp_package parsnp

# Or using assembly directory (requires sample name matching)
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe/output \
  --input /path/to/assembly/files \
  --outdir results \
  --snp_package parsnp
```

### Parameters

#### Required Parameters

- `--poppipe_output`: Path to PopPIPE-bp output directory containing:
  - Cluster assignments file (preferably `all_clusters.txt`, or fallback to `*_clusters.csv`)
  - Strain directories under `strain/` or `strains/` folder
  - Optional: backbone tree files for grafting

- `--input`: Path to TSV file with sample-to-file mapping (e.g., `rfiles.txt`, `combined.txt`) **OR** directory containing assembly files

- `--outdir`: Output directory for results

#### Optional Parameters

- `--ref`: Reference genome file (if not provided, largest genome in each cluster will be used)
- `--min_cluster_size`: Minimum number of samples required in a cluster for analysis (default: 3)
- `--enable_tree_grafting`: Enable/disable tree grafting (default: true)
- `--snp_package`: SNP calling method (default: "parsnp")

### Input Formats

#### Option 1: TSV File Input (Recommended)

Use the same TSV file that PopPIPE uses to identify samples. This file should contain tab-separated values with sample names in the first column and file paths in the second column:

```
sample1	/path/to/sample1.fasta
sample2	/path/to/sample2.fasta
sample3	/path/to/sample3.fasta
```

Common file names include:
- `rfiles.txt` (PopPUNK rfile format)
- `combined.txt` (combined sample list)
- Any tab-separated file with sample-name and file-path columns

#### Option 2: Directory Input

Provide a directory containing assembly files. The pipeline will attempt to match filenames to sample names in the PopPIPE cluster assignments.

### Expected PopPIPE-bp Output Structure

The pipeline expects PopPIPE-bp output with the following structure:

```
poppipe_output/
├── all_clusters.txt               # Cluster assignments (preferred format)
├── strain/                        # Strain-specific directories
│   ├── cluster_1/
│   │   ├── njtree.nwk            # NJ tree for scaling
│   │   ├── besttree.nwk          # ML tree
│   │   ├── fastbaps_clusters.txt
│   │   └── ...
│   ├── cluster_2/
│   │   └── ...
│   └── other/
│       └── njtree.nwk            # Backbone tree
└── full_tree.nwk                 # Optional: overall tree
```

## Workflow Steps

### 1. PopPIPE Output Parsing

The `PARSE_POPPIPE_OUTPUT` module:
- Identifies cluster assignment CSV files
- Maps sample names to assembly files
- Creates cluster-specific file lists
- Generates cluster summary information

### 2. Cluster-based SNP Analysis

For each cluster with ≥ `min_cluster_size` samples:
- Runs ParSNP for core genome alignment
- Calculates pairwise SNP distances
- Performs recombination detection (Gubbins)
- Builds phylogenetic trees
- Creates distance matrices

### 3. Tree Grafting

The `TREE_GRAFT` module implements the PopPIPE algorithm:
- Loads cluster-specific ML trees from SNP analysis
- Loads PopPIPE NJ trees for branch length scaling
- Scales ML tree branches based on NJ tree lengths
- Grafts scaled ML trees onto backbone tree
- Outputs final integrated phylogeny

## Output Files

### Main Outputs

- `grafted_tree.nwk`: Final integrated phylogenetic tree
- `cluster_summary.tsv`: Summary of clusters and sample counts
- `tree_graft_log.txt`: Log of tree grafting process

### Per-cluster Outputs

For each cluster, the pipeline generates:
- `cluster_X/`: Directory containing cluster-specific results
  - SNP distance matrices
  - Phylogenetic trees
  - Core genome alignments
  - Recombination analysis results

### Summary Files

- `Summary.QC_File_Checks.tsv`: Quality control checks
- `software_versions.yml`: Software version information

## Example Usage Scenarios

### Scenario 1: Using TSV File (Recommended)

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input poppipe_results/rfiles.txt \
  --outdir snp_analysis_results/
```

### Scenario 2: Using Directory Input

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input assemblies/ \
  --outdir snp_analysis_results/
```

### Scenario 3: With Reference Genome

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input poppipe_results/rfiles.txt \
  --ref reference_genome.fasta \
  --outdir snp_analysis_results/
```

### Scenario 4: Custom Cluster Size

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input poppipe_results/rfiles.txt \
  --min_cluster_size 5 \
  --outdir snp_analysis_results/
```

## Troubleshooting

### Common Issues

1. **Sample Name Matching**: If samples aren't matched between PopPIPE output and assembly files:
   - **Recommended solution**: Use TSV input (`--input rfiles.txt`) instead of directory input
   - TSV input eliminates sample name matching issues by using exact file paths
   - For directory input: Check that assembly filenames correspond to sample names in the CSV
   - The pipeline attempts fuzzy matching for common naming variations when using directory input

2. **Missing PopPIPE Data**: If PopPIPE strain directories are missing:
   - Tree grafting will create a combined tree from cluster analyses
   - Check the `tree_graft_log.txt` for details

3. **Small Clusters**: Clusters below `min_cluster_size` are skipped:
   - Adjust `--min_cluster_size` parameter if needed
   - Check `cluster_summary.tsv` for cluster sizes

### File Format Requirements

- **Assembly files**: FASTA format with extensions: `.fasta`, `.fas`, `.fna`, `.fsa`, `.fa` (with optional `.gz` compression)
- **PopPIPE CSV**: Must contain columns for sample names and cluster assignments
- **Tree files**: Newick format (`.nwk`, `.tree`, `.treefile`)

## Integration Benefits

This integration provides several advantages:

1. **Scalability**: Analyze large datasets by processing clusters separately
2. **Accuracy**: Detailed SNP analysis within homogeneous clusters
3. **Comprehensive Results**: Integrated phylogeny representing all samples
4. **PopPIPE Compatibility**: Leverages existing PopPIPE clustering and tree structure
5. **Flexible Input**: Works with various PopPIPE output configurations

## Citation

If you use this integration, please cite both:
- The original wf-assembly-snps pipeline
- PopPIPE-bp: McHugh, M. P., et al. (2025). Integrated population clustering and genomic epidemiology with PopPIPE. Microbial Genomics.