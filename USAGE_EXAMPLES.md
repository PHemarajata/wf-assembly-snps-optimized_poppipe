# PopPIPE-bp Integration Usage Examples

This document provides detailed examples of how to use the PopPIPE-bp integration with different input formats.

## Input Format Options

### Option 1: TSV File Input (Recommended)

The most reliable method is to use the same TSV file that PopPIPE uses for sample identification. This eliminates any sample name matching issues.

#### TSV File Format

The TSV file should contain tab-separated values with:
- **Column 1**: Sample name (exactly as it appears in the PopPIPE cluster CSV)
- **Column 2**: Full path to the assembly file

Example `rfiles.txt`:
```
sample001	/data/assemblies/sample001.fasta
sample002	/data/assemblies/sample002.fasta
sample003	/data/assemblies/sample003.fasta
GCA_123456.1	/data/assemblies/GCA_123456.1_genomic.fasta
GCA_789012.1	/data/assemblies/GCA_789012.1_genomic.fasta
```

#### Usage with TSV File

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe/output \
  --input /path/to/poppipe/rfiles.txt \
  --outdir results_tsv \
  --snp_package parsnp
```

### Option 2: Directory Input

If you don't have a TSV file, you can provide a directory containing assembly files. The pipeline will attempt to match filenames to sample names.

#### Directory Structure

```
assemblies/
├── sample001.fasta
├── sample002.fasta
├── sample003.fasta
├── GCA_123456.1_genomic.fasta
└── GCA_789012.1_genomic.fasta
```

#### Usage with Directory

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output /path/to/poppipe/output \
  --input /path/to/assemblies \
  --outdir results_dir \
  --snp_package parsnp
```

## Complete Examples

### Example 1: Basic PopPIPE Integration with TSV

```bash
# Assuming you have PopPIPE output and the original rfiles.txt
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input poppipe_results/rfiles.txt \
  --outdir snp_analysis_results/ \
  --min_cluster_size 3
```

**Expected PopPIPE structure:**
```
poppipe_results/
├── combined_clusters.csv
├── rfiles.txt
├── strain/
│   ├── cluster_1/
│   ├── cluster_2/
│   └── other/
└── full_tree.nwk
```

### Example 2: With Reference Genome

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input poppipe_results/rfiles.txt \
  --ref reference_genome.fasta \
  --outdir snp_analysis_with_ref/ \
  --min_cluster_size 5
```

### Example 3: Using Combined TSV File

If you have a different TSV file (e.g., `combined.txt`):

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output poppipe_results/ \
  --input poppipe_results/combined.txt \
  --outdir snp_analysis_combined/ \
  --enable_tree_grafting true
```

### Example 4: Directory Input with Custom Parameters

```bash
nextflow run main.nf \
  -profile singularity \
  --poppipe_output poppipe_results/ \
  --input /data/assemblies/ \
  --outdir snp_analysis_dir/ \
  --min_cluster_size 4 \
  --create_excel_outputs true
```

### Example 5: Large Dataset with Resource Optimization

```bash
nextflow run main.nf \
  -profile docker \
  --poppipe_output large_poppipe_results/ \
  --input large_poppipe_results/rfiles.txt \
  --outdir large_snp_analysis/ \
  --min_cluster_size 10 \
  --max_cpus 32 \
  --max_memory 200.GB
```

## Creating TSV Files

### From PopPUNK rfiles

If you have the original PopPUNK rfiles.txt, you can use it directly:

```bash
# PopPUNK rfiles.txt format is already compatible
cp poppunk_rfiles.txt poppipe_input.txt
```

### From Directory Listing

If you need to create a TSV file from a directory:

```bash
# Create TSV file from directory
cd /path/to/assemblies
for file in *.fasta; do
    sample_name=$(basename "$file" .fasta)
    echo -e "${sample_name}\t$(pwd)/${file}"
done > rfiles.txt
```

### From Samplesheet

If you have a CSV samplesheet:

```bash
# Convert CSV samplesheet to TSV format
tail -n +2 samplesheet.csv | \
awk -F',' '{print $1 "\t" $2}' > rfiles.txt
```

## Troubleshooting Input Issues

### TSV File Issues

1. **File paths not found**:
   ```bash
   # Check if paths in TSV are absolute and correct
   while IFS=$'\t' read -r sample path; do
       if [[ ! -f "$path" ]]; then
           echo "Missing: $sample -> $path"
       fi
   done < rfiles.txt
   ```

2. **Sample name mismatches**:
   ```bash
   # Compare sample names in TSV vs PopPIPE CSV
   cut -f1 rfiles.txt | sort > tsv_samples.txt
   cut -d',' -f1 combined_clusters.csv | tail -n +2 | sort > csv_samples.txt
   diff tsv_samples.txt csv_samples.txt
   ```

### Directory Input Issues

1. **Sample name matching problems**:
   - Use TSV input instead for exact matching
   - Check that filenames correspond to sample names in PopPIPE CSV
   - Remove common suffixes like `_contigs`, `_genomic`, `-SPAdes`

2. **File extension issues**:
   - Ensure files have supported extensions: `.fasta`, `.fas`, `.fna`, `.fsa`, `.fa`
   - Compressed files (`.gz`) are supported

## Output Interpretation

### Successful Run Outputs

```
results/
├── grafted_tree.nwk              # Final integrated phylogeny
├── cluster_summary.tsv           # Cluster information
├── tree_graft_log.txt           # Grafting process log
├── cluster_1/                   # Results for cluster 1
│   ├── distance_matrix.tsv
│   ├── phylogenetic_tree.nwk
│   └── core_alignment.fasta
├── cluster_2/                   # Results for cluster 2
│   └── ...
└── Summaries/
    └── Summary.QC_File_Checks.tsv
```

### Key Files to Check

1. **`cluster_summary.tsv`**: Verify all expected clusters are present with correct sample counts
2. **`tree_graft_log.txt`**: Check for any grafting issues or warnings
3. **`Summary.QC_File_Checks.tsv`**: Review quality control results

## Performance Considerations

### For Large Datasets

- Use TSV input for better performance (no filename matching required)
- Increase `--min_cluster_size` to focus on larger clusters
- Use appropriate resource limits (`--max_cpus`, `--max_memory`)
- Consider using Singularity for better resource management

### For Small Datasets

- Default parameters usually work well
- Directory input is acceptable for small numbers of samples
- Consider lowering `--min_cluster_size` to include smaller clusters

## Integration with PopPIPE Workflow

### Typical Workflow

1. **Run PopPUNK** to generate initial clusters
2. **Run PopPIPE** to refine clusters and generate strain-specific analyses
3. **Run this integration** to perform detailed SNP analysis per cluster
4. **Use grafted tree** for comprehensive phylogenetic analysis

### File Reuse

- Use the same `rfiles.txt` from PopPUNK/PopPIPE
- Reuse PopPIPE output directory directly
- Reference genomes can be the same as used in PopPIPE

This integration seamlessly extends the PopPIPE workflow with detailed SNP analysis while maintaining compatibility with existing PopPIPE outputs and file formats.