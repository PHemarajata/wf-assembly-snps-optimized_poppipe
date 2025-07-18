# PopPUNK Cluster-Based SNP Analysis

This pipeline has been enhanced to support PopPUNK cluster-based processing, allowing you to run SNP analysis on each PopPUNK cluster separately instead of processing all samples together. This is particularly useful when dealing with large datasets where computational resources are limited.

## Overview

When PopPUNK clustering is enabled, the pipeline will:

1. Parse the PopPUNK cluster assignments
2. Group samples by cluster
3. Run SNP analysis independently on each cluster
4. Generate separate outputs for each cluster
5. Combine summary statistics across all clusters

## Requirements

### PopPUNK Results File

You need a CSV file with PopPUNK cluster assignments containing at least these columns:
- `Taxon`: Sample name (should match your FASTA file names)
- `Cluster`: Cluster ID/number

Example format:
```csv
Taxon,Cluster
sample001,1
sample002,1
sample003,2
sample004,2
sample005,3
```

### Input Files

Your FASTA files should be named consistently with the sample names in the PopPUNK results. The pipeline supports these extensions:
- `.fasta`, `.fas`, `.fna`, `.fsa`, `.fa`
- Gzipped versions: `.fasta.gz`, `.fas.gz`, etc.

## Usage

### Basic Command

```bash
nextflow run main.nf \
  -profile gcp_vm \
  --input /path/to/fasta/files \
  --poppunk_clusters /path/to/full_assign_final.csv \
  --outdir results_clustered \
  --snp_package parsnp
```

### Parameters

- `--poppunk_clusters`: Path to PopPUNK cluster assignment CSV file
- `--min_cluster_size`: Minimum number of samples required per cluster (default: 3)
- `--input`: Directory containing FASTA files
- `--outdir`: Output directory for results

### Preparing Your Data

Use the provided helper script to validate and prepare your PopPUNK data:

```bash
python bin/prepare_poppunk_data.py \
  --poppunk_csv /path/to/full_assign_final.csv \
  --input_dir /path/to/fasta/files \
  --output_dir poppunk_prepared \
  --min_cluster_size 3
```

This script will:
- Validate your PopPUNK results file
- Check which samples have corresponding FASTA files
- Filter clusters by minimum size
- Generate a summary report
- Create an example command to run the pipeline

## Output Structure

When running in cluster mode, outputs are organized by cluster:

```
results_clustered/
├── Parsnp/
│   ├── cluster_1/
│   │   ├── Parsnp.ggr
│   │   ├── Parsnp.tree
│   │   └── ...
│   ├── cluster_2/
│   │   ├── Parsnp.ggr
│   │   ├── Parsnp.tree
│   │   └── ...
│   └── ...
├── Summaries/
│   ├── Summary.Parsnp.Distance_Matrix.cluster_1.tsv
│   ├── Summary.Parsnp.Distance_Matrix.cluster_2.tsv
│   └── ...
└── pipeline_info/
    └── ...
```

## Performance Considerations

### GCP VM Optimization

The pipeline has been optimized for a 64-CPU, 412GB RAM GCP VM:

- **Resource allocation**: Up to 60 CPUs and 400GB RAM per process
- **Parallel processing**: Multiple clusters can be processed simultaneously
- **Memory management**: Efficient memory usage for large datasets

### Cluster Size Recommendations

- **Minimum cluster size**: 3-5 samples (configurable)
- **Optimal cluster size**: 10-50 samples per cluster
- **Large clusters**: Clusters with >100 samples may require more resources

## Troubleshooting

### Common Issues

1. **Sample name mismatch**: Ensure sample names in PopPUNK CSV match FASTA filenames
2. **Small clusters**: Adjust `--min_cluster_size` parameter
3. **Memory issues**: Reduce the number of concurrent processes or increase VM resources

### Validation

Before running the full pipeline, use the preparation script to validate your data:

```bash
python bin/prepare_poppunk_data.py --poppunk_csv your_file.csv --input_dir your_fasta_dir
```

### Monitoring

Monitor resource usage during execution:
- Check CPU utilization: `htop`
- Check memory usage: `free -h`
- Monitor disk space: `df -h`

## Example Workflow

1. **Prepare data**:
   ```bash
   python bin/prepare_poppunk_data.py \
     --poppunk_csv poppunk_results/full_assign_final.csv \
     --input_dir genome_assemblies/ \
     --min_cluster_size 5
   ```

2. **Run pipeline**:
   ```bash
   nextflow run main.nf \
     -profile gcp_vm \
     --input genome_assemblies/ \
     --poppunk_clusters poppunk_prepared/filtered_poppunk_clusters.csv \
     --min_cluster_size 5 \
     --outdir results_by_cluster
   ```

3. **Check results**:
   ```bash
   ls results_by_cluster/Summaries/
   ```

## Performance Benchmarks

On a 64-CPU, 412GB RAM GCP VM:
- **Small clusters** (3-10 samples): ~5-15 minutes per cluster
- **Medium clusters** (10-50 samples): ~15-60 minutes per cluster  
- **Large clusters** (50+ samples): ~1-4 hours per cluster

Multiple clusters can be processed in parallel, significantly reducing total runtime compared to processing all samples together.