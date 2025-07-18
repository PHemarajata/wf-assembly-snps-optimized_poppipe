#!/usr/bin/env python3
"""
Analyze PopPUNK cluster data and provide recommendations for SNP analysis pipeline.
"""

import pandas as pd
from pathlib import Path

def analyze_poppunk_clusters(csv_file):
    """Analyze PopPUNK cluster assignments and provide detailed statistics."""
    
    # Read the data
    df = pd.read_csv(csv_file)
    
    print("=" * 60)
    print("PopPUNK Cluster Analysis Report")
    print("=" * 60)
    print(f"Input file: {csv_file}")
    print(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Basic statistics
    total_samples = len(df)
    total_clusters = df['Cluster'].nunique()
    
    print("BASIC STATISTICS:")
    print(f"  Total samples: {total_samples}")
    print(f"  Total clusters: {total_clusters}")
    print(f"  Average samples per cluster: {total_samples/total_clusters:.1f}")
    print()
    
    # Cluster size distribution
    cluster_counts = df['Cluster'].value_counts().sort_index()
    
    print("CLUSTER SIZE DISTRIBUTION:")
    size_bins = [1, 2, 3, 5, 10, 20, 50, 100, float('inf')]
    size_labels = ['1', '2', '3-4', '5-9', '10-19', '20-49', '50-99', '100+']
    
    for i, (lower, upper) in enumerate(zip(size_bins[:-1], size_bins[1:])):
        if upper == float('inf'):
            count = sum(cluster_counts >= lower)
            samples = sum(cluster_counts[cluster_counts >= lower])
        else:
            count = sum((cluster_counts >= lower) & (cluster_counts < upper))
            samples = sum(cluster_counts[(cluster_counts >= lower) & (cluster_counts < upper)])
        
        if count > 0:
            print(f"  {size_labels[i]} samples: {count} clusters ({samples} total samples)")
    print()
    
    # Recommendations for different minimum cluster sizes
    print("PROCESSING RECOMMENDATIONS:")
    for min_size in [3, 5, 10, 15, 20]:
        eligible_clusters = cluster_counts[cluster_counts >= min_size]
        if len(eligible_clusters) > 0:
            total_samples_eligible = eligible_clusters.sum()
            largest_cluster = eligible_clusters.max()
            print(f"  Min cluster size {min_size}: {len(eligible_clusters)} clusters, {total_samples_eligible} samples")
            print(f"    Largest cluster: {largest_cluster} samples")
            print(f"    Coverage: {total_samples_eligible/total_samples*100:.1f}% of all samples")
        else:
            print(f"  Min cluster size {min_size}: No eligible clusters")
        print()
    
    # Detailed cluster breakdown
    print("DETAILED CLUSTER BREAKDOWN:")
    print("Cluster ID | Sample Count | Percentage")
    print("-" * 40)
    for cluster_id, count in cluster_counts.items():
        percentage = count / total_samples * 100
        print(f"    {cluster_id:2d}     |     {count:3d}      |   {percentage:5.1f}%")
    print()
    
    # Sample names analysis
    print("SAMPLE NAME PATTERNS:")
    sample_prefixes = {}
    for sample in df['Taxon']:
        # Extract prefix patterns
        if '_' in sample:
            prefix = sample.split('_')[0]
        elif '-' in sample:
            prefix = sample.split('-')[0]
        else:
            prefix = sample[:10] if len(sample) > 10 else sample
        
        sample_prefixes[prefix] = sample_prefixes.get(prefix, 0) + 1
    
    # Show most common prefixes
    common_prefixes = sorted(sample_prefixes.items(), key=lambda x: x[1], reverse=True)[:10]
    print("  Most common sample prefixes:")
    for prefix, count in common_prefixes:
        print(f"    {prefix}: {count} samples")
    print()
    
    # Computational estimates
    print("COMPUTATIONAL ESTIMATES (64-CPU, 412GB RAM VM):")
    print()
    
    # Regular processing estimate
    print("  Regular processing (all samples together):")
    if total_samples <= 50:
        time_est = "2-6 hours"
    elif total_samples <= 200:
        time_est = "6-24 hours"
    else:
        time_est = "24-72 hours"
    print(f"    Estimated time: {time_est}")
    print(f"    Memory usage: High (may require process_max resources)")
    print()
    
    # Cluster processing estimates
    print("  Cluster-based processing:")
    for min_size in [3, 5, 10]:
        eligible_clusters = cluster_counts[cluster_counts >= min_size]
        if len(eligible_clusters) > 0:
            # Estimate time per cluster
            cluster_times = []
            for size in eligible_clusters:
                if size <= 10:
                    cluster_times.append(15)  # 15 minutes
                elif size <= 50:
                    cluster_times.append(60)  # 1 hour
                else:
                    cluster_times.append(240)  # 4 hours
            
            # Assuming parallel processing
            max_parallel = min(8, len(eligible_clusters))  # Conservative estimate
            total_time_parallel = max(cluster_times) + sum(sorted(cluster_times, reverse=True)[max_parallel:]) / max_parallel
            
            print(f"    Min cluster size {min_size}:")
            print(f"      Clusters to process: {len(eligible_clusters)}")
            print(f"      Estimated total time: {total_time_parallel/60:.1f} hours (parallel)")
            print(f"      Memory per cluster: Moderate to high")
            print(f"      Recommended for: {'Yes' if len(eligible_clusters) >= 3 else 'No'}")
            print()
    
    # Generate cluster assignment files for different minimum sizes
    output_dir = Path("cluster_analysis_output")
    output_dir.mkdir(exist_ok=True)
    
    for min_size in [3, 5, 10]:
        eligible_clusters = cluster_counts[cluster_counts >= min_size]
        if len(eligible_clusters) > 0:
            filtered_df = df[df['Cluster'].isin(eligible_clusters.index)]
            output_file = output_dir / f"poppunk_clusters_min{min_size}.csv"
            filtered_df.to_csv(output_file, index=False)
            print(f"Saved filtered clusters (min size {min_size}) to: {output_file}")
    
    # Generate summary statistics file
    summary_file = output_dir / "cluster_analysis_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("PopPUNK Cluster Analysis Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Total samples: {total_samples}\n")
        f.write(f"Total clusters: {total_clusters}\n")
        f.write(f"Average cluster size: {total_samples/total_clusters:.1f}\n\n")
        
        f.write("Cluster size distribution:\n")
        for cluster_id, count in cluster_counts.items():
            f.write(f"  Cluster {cluster_id}: {count} samples\n")
    
    print(f"Saved detailed summary to: {summary_file}")
    print()
    
    # Recommendations
    print("FINAL RECOMMENDATIONS:")
    print()
    
    # Find optimal minimum cluster size
    best_min_size = None
    best_coverage = 0
    for min_size in [3, 5, 10, 15]:
        eligible_clusters = cluster_counts[cluster_counts >= min_size]
        if len(eligible_clusters) >= 3:  # At least 3 clusters
            coverage = eligible_clusters.sum() / total_samples
            if coverage > best_coverage:
                best_coverage = coverage
                best_min_size = min_size
    
    if best_min_size:
        print(f"  RECOMMENDED: Use minimum cluster size of {best_min_size}")
        eligible_clusters = cluster_counts[cluster_counts >= best_min_size]
        print(f"    - Processes {len(eligible_clusters)} clusters")
        print(f"    - Covers {eligible_clusters.sum()} samples ({best_coverage*100:.1f}% of total)")
        print(f"    - Estimated runtime: Much faster than regular processing")
        print()
        
        print("  COMMAND TO RUN:")
        print(f"    nextflow run main.nf \\")
        print(f"      -profile gcp_vm \\")
        print(f"      --input /path/to/your/fasta/files \\")
        print(f"      --poppunk_clusters cluster_analysis_output/poppunk_clusters_min{best_min_size}.csv \\")
        print(f"      --min_cluster_size {best_min_size} \\")
        print(f"      --outdir results_clustered \\")
        print(f"      --snp_package parsnp")
    else:
        print("  RECOMMENDED: Use regular processing (clusters too small for efficient parallel processing)")
        print("    nextflow run main.nf \\")
        print("      -profile gcp_vm \\")
        print("      --input /path/to/your/fasta/files \\")
        print("      --outdir results_regular \\")
        print("      --snp_package parsnp")
    
    print()
    print("=" * 60)
    print("Analysis complete!")
    print("=" * 60)

if __name__ == "__main__":
    csv_file = "results_poppunk/poppunk_full/full_assign_final.csv"
    analyze_poppunk_clusters(csv_file)