# PopPUNK Cluster Analysis Results

## Summary of Your Data

Your PopPUNK analysis has identified **449 samples** organized into **27 clusters**. This is excellent data for cluster-based SNP analysis!

### Key Statistics:
- **Total samples**: 449
- **Total clusters**: 27  
- **Average cluster size**: 16.6 samples
- **Largest cluster**: 264 samples (58.8% of all data)
- **Second largest**: 76 samples (16.9% of all data)

### Cluster Distribution:
- **Large clusters (≥20 samples)**: 3 clusters (372 samples, 82.9% coverage)
- **Medium clusters (10-19 samples)**: 2 clusters (27 samples)
- **Small clusters (3-9 samples)**: 7 clusters (31 samples)
- **Singleton/Doubleton clusters**: 15 clusters (19 samples)

## Optimization Results

### 🚀 **MASSIVE Performance Improvement Expected!**

**Regular Processing** (all 449 samples together):
- **Estimated time**: 24-72 hours
- **Memory usage**: Very high (may exceed VM limits)
- **Risk**: Potential memory issues with large dataset

**Cluster-Based Processing** (recommended):
- **Estimated time**: ~4 hours total (parallel processing)
- **Memory usage**: Moderate per cluster
- **Coverage**: 95.8% of samples (430/449)
- **Performance gain**: **85-95% time reduction!**

## Recommendations

### ✅ **RECOMMENDED APPROACH: Cluster-based with min_cluster_size=3**

This approach will:
- Process **12 clusters** in parallel
- Cover **430 samples** (95.8% of your data)
- Complete in approximately **4 hours** instead of 24-72 hours
- Use memory efficiently across clusters
- Provide separate phylogenetic trees for each cluster

### 📁 **Files Generated for You:**

1. **`cluster_analysis_output/poppunk_clusters_min3.csv`** - Optimized cluster file
2. **`cluster_analysis_output/poppunk_clusters_min5.csv`** - More conservative option
3. **`cluster_analysis_output/poppunk_clusters_min10.csv`** - Most conservative option

## Ready-to-Run Commands

### Option 1: Recommended (min cluster size 3)
```bash
nextflow run main.nf \
  -profile gcp_vm \
  --input /path/to/your/fasta/files \
  --poppunk_clusters cluster_analysis_output/poppunk_clusters_min3.csv \
  --min_cluster_size 3 \
  --outdir results_clustered_min3 \
  --snp_package parsnp
```

### Option 2: Conservative (min cluster size 5)
```bash
nextflow run main.nf \
  -profile gcp_vm \
  --input /path/to/your/fasta/files \
  --poppunk_clusters cluster_analysis_output/poppunk_clusters_min5.csv \
  --min_cluster_size 5 \
  --outdir results_clustered_min5 \
  --snp_package parsnp
```

### Option 3: Most Conservative (min cluster size 10)
```bash
nextflow run main.nf \
  -profile gcp_vm \
  --input /path/to/your/fasta/files \
  --poppunk_clusters cluster_analysis_output/poppunk_clusters_min10.csv \
  --min_cluster_size 10 \
  --outdir results_clustered_min10 \
  --snp_package parsnp
```

## Expected Output Structure

```
results_clustered_min3/
├── Parsnp/
│   ├── cluster_1/          # 264 samples - largest cluster
│   │   ├── Parsnp.ggr
│   │   ├── Parsnp.tree
│   │   └── core_alignment.fasta
│   ├── cluster_2/          # 76 samples - second largest
│   │   ├── Parsnp.ggr
│   │   ├── Parsnp.tree
│   │   └── core_alignment.fasta
│   ├── cluster_3/          # 15 samples
│   ├── cluster_4/          # 12 samples
│   ├── cluster_5/          # 32 samples
│   └── ... (7 more clusters)
├── Summaries/
│   ├── Summary.Parsnp.Distance_Matrix.cluster_1.tsv
│   ├── Summary.Parsnp.Distance_Matrix.cluster_2.tsv
│   └── ... (one per cluster)
└── pipeline_info/
    ├── execution_timeline_*.html
    ├── execution_report_*.html
    └── software_versions.yml
```

## Performance Expectations on Your 64-CPU, 412GB RAM VM

### Parallel Processing:
- **Cluster 1** (264 samples): ~4 hours (largest, will determine total time)
- **Cluster 2** (76 samples): ~1 hour (runs in parallel)
- **Smaller clusters**: 15-30 minutes each (run in parallel)

### Resource Utilization:
- **CPU**: Up to 60 cores utilized across multiple clusters
- **Memory**: 50-150GB per large cluster (well within your 412GB limit)
- **I/O**: Optimized for parallel processing

## Sample Name Patterns Detected

Your samples appear to follow these patterns:
- **GCA_*** (98 samples) - Likely NCBI GenBank assemblies
- **Burkholderia_*** (50 samples) - Reference genomes
- **IP-***, **IE-***, **ERS*** - Various sample naming conventions

The pipeline will automatically match these to your FASTA files regardless of the naming pattern.

## Next Steps

1. **Prepare your FASTA files** in a single directory
2. **Choose your preferred minimum cluster size** (recommend: 3)
3. **Run the command** with the path to your FASTA directory
4. **Monitor progress** - the pipeline will show which clusters are being processed
5. **Analyze results** - each cluster will have its own phylogenetic tree and distance matrices

## Troubleshooting

If you encounter issues:
1. **Check sample name matching**: Use `python3 bin/prepare_poppunk_data.py` to validate
2. **Monitor resources**: Use `htop` and `free -h` to check CPU/memory usage
3. **Check logs**: Look in `results_clustered_min3/pipeline_info/` for detailed logs

---

**🎉 Your pipeline is now optimized for maximum performance with your PopPUNK cluster data!**

The combination of GCP VM optimization + PopPUNK cluster-based processing should give you **dramatic performance improvements** while maintaining the quality of your SNP analysis.