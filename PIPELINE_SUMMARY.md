# Modified wf-assembly-snps Pipeline with PopPIPE-bp Integration

## Overview
This is a modified version of the wf-assembly-snps pipeline that includes integration with PopPIPE-bp outputs for cluster-based SNP analysis and tree grafting.

## Key Modifications

### New Features Added
1. **PopPIPE-bp Integration**: Parse PopPIPE-bp cluster outputs and perform SNP analysis per cluster
2. **Tree Grafting**: Combine cluster-specific phylogenies using the PopPIPE algorithm
3. **Cluster-based Analysis**: Process large datasets by analyzing clusters separately

### New Files Added
- `modules/local/parse_poppipe_output/main.nf` - PopPIPE output parser
- `modules/local/tree_graft/main.nf` - Tree grafting implementation
- `workflows/poppipe_snp_analysis.nf` - PopPIPE integration workflow
- `POPPIPE_INTEGRATION.md` - Comprehensive usage guide
- `test_poppipe_integration.nf` - Test script
- `MODIFICATION_SUMMARY.md` - Detailed change summary

### Modified Files
- `main.nf` - Added PopPIPE workflow option
- `conf/params.config` - Added PopPIPE parameters
- `nextflow_schema.json` - Updated parameter schema
- `README.md` - Added PopPIPE usage examples

## Usage

### Standard SNP Analysis
```bash
nextflow run main.nf \
  -profile docker \
  --input assemblies/ \
  --outdir results/
```

### PopPIPE-bp Integration
```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe/results \
  --input /path/to/assemblies \
  --outdir results/
```

### PopPUNK Cluster Analysis
```bash
nextflow run main.nf \
  -profile docker \
  --poppunk_clusters clusters.csv \
  --input assemblies/ \
  --outdir results/
```

## Key Parameters

### PopPIPE Integration
- `--poppipe_output`: Path to PopPIPE-bp output directory
- `--min_cluster_size`: Minimum samples per cluster (default: 3)
- `--enable_tree_grafting`: Enable tree grafting (default: true)

### Standard Parameters
- `--input`: Input assembly directory or samplesheet
- `--ref`: Reference genome (optional)
- `--outdir`: Output directory
- `--snp_package`: SNP calling method (default: parsnp)

## Pipeline Size
- Total size: ~4.2MB (cleaned)
- Nextflow files: 26 files
- Ready for efficient download and deployment

## Requirements
- Nextflow >=22.04.3
- Docker or Singularity
- 8GB+ RAM recommended
- 8+ CPUs recommended

## Support
See `POPPIPE_INTEGRATION.md` for detailed usage instructions and troubleshooting.