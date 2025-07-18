# PopPIPE-bp all_clusters.txt Integration Guide

This document provides comprehensive guidance on using the PopPIPE-bp integration with the `all_clusters.txt` format.

## Overview

The pipeline now fully supports PopPIPE-bp's `all_clusters.txt` output format, which contains detailed cluster and subcluster information with the following structure:

```csv
Taxon,Strain,Subcluster_1,Subcluster_2,Subcluster_3
sample001,1,1_1,1_1_1,1_1_1_1
sample002,1,1_1,1_1_1,1_1_1_2
sample003,1,1_2,1_2_1,1_2_1_1
sample004,2,2_1,2_1_1,2_1_1_1
sample005,2,2_1,2_1_2,2_1_2_1
sample006,2,2_2,2_2_1,2_2_1_1
```

## Column Definitions

- **Taxon**: Sample identifier (exactly as used in assembly files)
- **Strain**: Primary cluster/strain assignment (used for SNP analysis grouping)
- **Subcluster_1**: First level of subclustering within strain
- **Subcluster_2**: Second level of subclustering
- **Subcluster_3**: Third level of subclustering

## Expected PopPIPE-bp Output Structure

The pipeline expects your PopPIPE-bp output to have this structure:

```
poppipe_output/
├── all_clusters.txt                    # Main cluster assignment file
├── strains/                           # Strain-specific directories
│   ├── 1/                            # Strain 1 directory
│   │   ├── align_variants.aln        # Variant alignment
│   │   ├── besttree.nwk             # ML phylogenetic tree
│   │   ├── dists.npy                # Distance matrix (numpy)
│   │   ├── dists.pkl                # Distance matrix (pickle)
│   │   ├── fastbaps_clusters.txt    # FastBAPS subclusters
│   │   ├── names.txt                # Sample names
│   │   ├── njtree.nwk              # Neighbor-joining tree
│   │   ├── rfile.txt               # Sample-to-file mapping
│   │   └── split_kmers.skf         # SKA k-mer file
│   ├── 2/                          # Strain 2 directory
│   │   └── ... (same files as strain 1)
│   └── other/                      # Backbone/overall tree
│       └── njtree.nwk             # Overall NJ tree
└── rfiles.txt                     # Optional: sample-to-file mapping
```

## Usage Examples

### Basic Usage with all_clusters.txt

```bash
# Using all_clusters.txt with TSV file for assembly paths
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe/output \
  --input /path/to/poppipe/rfiles.txt \
  --outdir results_all_clusters \
  --snp_package parsnp
```

### Using Directory Input

```bash
# Using all_clusters.txt with assembly directory
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe/output \
  --input /path/to/assemblies \
  --outdir results_all_clusters \
  --snp_package parsnp
```

### With Custom Parameters

```bash
# Advanced usage with custom cluster size and reference
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe/output \
  --input /path/to/poppipe/rfiles.txt \
  --ref reference_genome.fasta \
  --min_cluster_size 5 \
  --outdir results_custom \
  --enable_tree_grafting true
```

## File Priority and Detection

The pipeline searches for cluster files in this priority order:

1. **`all_clusters.txt`** (highest priority)
   - `poppipe_output/all_clusters.txt`
   - `poppipe_output/output/all_clusters.txt`

2. **CSV files** (medium priority)
   - `poppipe_output/*_clusters.csv`
   - `poppipe_output/*.csv`
   - `poppipe_output/output/*_clusters.csv`

3. **Other TXT files** (lowest priority)
   - `poppipe_output/*_clusters.txt`
   - `poppipe_output/output/*_clusters.txt`

## What the Pipeline Does

### 1. Cluster File Parsing
- Automatically detects and reads `all_clusters.txt`
- Identifies `Taxon` column for sample names
- Uses `Strain` column for primary clustering
- Preserves subcluster information for reference

### 2. Sample-to-File Mapping
- Maps samples from `all_clusters.txt` to assembly files
- Supports both TSV input (rfiles.txt) and directory scanning
- Performs exact matching for TSV input
- Uses fuzzy matching for directory input

### 3. Strain-Based SNP Analysis
- Groups samples by `Strain` column
- Runs ParSNP analysis on each strain separately
- Generates strain-specific phylogenetic trees
- Calculates SNP distance matrices per strain

### 4. Tree Grafting
- Loads PopPIPE NJ trees from strain directories
- Scales ML tree branch lengths based on NJ trees
- Grafts detailed strain trees onto backbone structure
- Produces integrated phylogeny with accurate SNP representation

## Output Files

### Main Results
```
results/
├── grafted_tree.nwk              # Final integrated phylogenetic tree
├── cluster_summary.tsv           # Summary of strains and sample counts
├── tree_graft_log.txt           # Tree grafting process log
└── cluster_info.json           # PopPIPE strain directory information
```

### Per-Strain Results
```
results/
├── cluster_1/                   # Results for strain 1
│   ├── distance_matrix.tsv      # SNP distance matrix
│   ├── masked_distance_matrix.tsv # Recombination-masked distances
│   ├── phylogenetic_tree.nwk    # ML phylogenetic tree
│   ├── core_alignment.fasta     # Core genome alignment
│   └── masked_alignment.fasta   # Recombination-masked alignment
├── cluster_2/                   # Results for strain 2
│   └── ... (same structure)
└── Summaries/
    └── Summary.QC_File_Checks.tsv # Quality control summary
```

## Validation and Testing

The integration has been tested with:
- ✅ Correct parsing of `all_clusters.txt` format
- ✅ Proper identification of Taxon and Strain columns
- ✅ Successful sample-to-file mapping
- ✅ Strain-based clustering (strain 1: 3 samples, strain 2: 3 samples)
- ✅ PopPIPE strain directory detection
- ✅ Subcluster information preservation

## Troubleshooting

### Common Issues

1. **all_clusters.txt not found**
   ```
   Error: No cluster file found in PopPIPE output directory
   ```
   - Ensure `all_clusters.txt` exists in the PopPIPE output directory
   - Check that the path to `--poppipe_output` is correct

2. **Sample matching failures**
   ```
   Warning: No match found for sample: sample_name
   ```
   - Use TSV input (rfiles.txt) for exact matching
   - Verify sample names in `all_clusters.txt` match assembly filenames

3. **Missing strain directories**
   ```
   Found 0 cluster directories in strain folder
   ```
   - Check that `strains/` directory exists in PopPIPE output
   - Verify strain directories contain expected files (njtree.nwk, besttree.nwk, etc.)

### Validation Commands

```bash
# Check all_clusters.txt format
head -5 /path/to/poppipe/output/all_clusters.txt

# Verify strain directories
ls -la /path/to/poppipe/output/strains/

# Check sample-to-file mapping
head -5 /path/to/poppipe/rfiles.txt
```

## Integration Benefits

1. **Hierarchical Clustering**: Supports PopPIPE's multi-level clustering (Strain → Subcluster_1 → Subcluster_2 → Subcluster_3)
2. **Accurate SNP Analysis**: Detailed SNP analysis within homogeneous strains
3. **Tree Grafting**: Combines strain-specific phylogenies using PopPIPE algorithm
4. **Scalability**: Processes large datasets by analyzing strains separately
5. **Compatibility**: Works seamlessly with existing PopPIPE-bp workflows

## Performance Recommendations

### For Large Datasets
- Use TSV input (rfiles.txt) for better performance
- Increase `--min_cluster_size` to focus on larger strains
- Use appropriate resource limits (`--max_cpus`, `--max_memory`)

### For Small Datasets
- Default parameters work well
- Consider lowering `--min_cluster_size` to include smaller strains
- Directory input is acceptable for small sample numbers

This integration provides a seamless bridge between PopPIPE-bp's population structure analysis and detailed SNP-based phylogenetic reconstruction.