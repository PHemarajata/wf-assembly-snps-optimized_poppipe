#!/usr/bin/env python3
"""
Script to validate and prepare PopPUNK cluster data for the SNP analysis pipeline.
"""

import pandas as pd
import argparse
import sys
from pathlib import Path
import os

def main():
    parser = argparse.ArgumentParser(description='Prepare PopPUNK cluster data for SNP analysis')
    parser.add_argument('--poppunk_csv', required=True, help='PopPUNK cluster assignment CSV file')
    parser.add_argument('--input_dir', required=True, help='Directory containing input FASTA files')
    parser.add_argument('--output_dir', default='poppunk_prepared', help='Output directory for prepared data')
    parser.add_argument('--min_cluster_size', type=int, default=3, help='Minimum cluster size to process')
    
    args = parser.parse_args()
    
    # Read PopPUNK results
    try:
        df = pd.read_csv(args.poppunk_csv)
        print(f"Loaded PopPUNK results with {len(df)} samples")
    except Exception as e:
        print(f"Error reading PopPUNK CSV file: {e}")
        sys.exit(1)
    
    # Check required columns
    required_cols = ['Taxon', 'Cluster']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)
    
    # Get input files
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"Error: Input directory does not exist: {input_dir}")
        sys.exit(1)
    
    # Find FASTA files
    fasta_extensions = ['.fasta', '.fas', '.fna', '.fsa', '.fa', '.fasta.gz', '.fas.gz', '.fna.gz', '.fsa.gz', '.fa.gz']
    fasta_files = []
    for ext in fasta_extensions:
        fasta_files.extend(list(input_dir.glob(f'*{ext}')))
    
    print(f"Found {len(fasta_files)} FASTA files in {input_dir}")
    
    # Create sample name to file mapping
    sample_to_file = {}
    for file_path in fasta_files:
        # Extract sample name from filename
        sample_name = file_path.stem
        # Remove common extensions
        for ext in ['.fasta', '.fas', '.fna', '.fsa', '.fa']:
            if sample_name.endswith(ext):
                sample_name = sample_name[:-len(ext)]
                break
        sample_to_file[sample_name] = file_path
    
    print(f"Mapped {len(sample_to_file)} sample names to files")
    
    # Check which samples from PopPUNK have corresponding files
    samples_with_files = []
    samples_without_files = []
    
    for _, row in df.iterrows():
        sample_name = row['Taxon']
        if sample_name in sample_to_file:
            samples_with_files.append(sample_name)
        else:
            samples_without_files.append(sample_name)
    
    print(f"Samples with files: {len(samples_with_files)}")
    print(f"Samples without files: {len(samples_without_files)}")
    
    if samples_without_files:
        print("Samples without corresponding files:")
        for sample in samples_without_files[:10]:  # Show first 10
            print(f"  - {sample}")
        if len(samples_without_files) > 10:
            print(f"  ... and {len(samples_without_files) - 10} more")
    
    # Filter dataframe to only include samples with files
    df_filtered = df[df['Taxon'].isin(samples_with_files)]
    
    # Analyze clusters
    cluster_stats = df_filtered.groupby('Cluster').size().sort_values(ascending=False)
    print(f"\nCluster statistics:")
    print(f"Total clusters: {len(cluster_stats)}")
    print(f"Clusters with >= {args.min_cluster_size} samples: {sum(cluster_stats >= args.min_cluster_size)}")
    
    print(f"\nTop 10 largest clusters:")
    for cluster, size in cluster_stats.head(10).items():
        print(f"  Cluster {cluster}: {size} samples")
    
    # Filter clusters by minimum size
    large_clusters = cluster_stats[cluster_stats >= args.min_cluster_size]
    df_final = df_filtered[df_filtered['Cluster'].isin(large_clusters.index)]
    
    print(f"\nFinal dataset:")
    print(f"Clusters to process: {len(large_clusters)}")
    print(f"Total samples to process: {len(df_final)}")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Save filtered results
    output_csv = output_dir / 'filtered_poppunk_clusters.csv'
    df_final.to_csv(output_csv, index=False)
    print(f"\nSaved filtered results to: {output_csv}")
    
    # Create summary report
    summary_file = output_dir / 'cluster_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("PopPUNK Cluster Analysis Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Input PopPUNK file: {args.poppunk_csv}\n")
        f.write(f"Input directory: {args.input_dir}\n")
        f.write(f"Total samples in PopPUNK: {len(df)}\n")
        f.write(f"Samples with FASTA files: {len(samples_with_files)}\n")
        f.write(f"Minimum cluster size: {args.min_cluster_size}\n")
        f.write(f"Clusters to process: {len(large_clusters)}\n")
        f.write(f"Total samples to process: {len(df_final)}\n\n")
        
        f.write("Cluster sizes:\n")
        for cluster, size in large_clusters.items():
            f.write(f"  Cluster {cluster}: {size} samples\n")
    
    print(f"Saved summary report to: {summary_file}")
    
    # Create example command
    example_cmd = f"""
# Example command to run the pipeline with PopPUNK clusters:
nextflow run main.nf \\
  -profile gcp_vm \\
  --input {args.input_dir} \\
  --poppunk_clusters {output_csv} \\
  --min_cluster_size {args.min_cluster_size} \\
  --outdir results_clustered \\
  --snp_package parsnp
"""
    
    print(example_cmd)
    
    with open(output_dir / 'example_command.sh', 'w') as f:
        f.write(example_cmd.strip())

if __name__ == '__main__':
    main()