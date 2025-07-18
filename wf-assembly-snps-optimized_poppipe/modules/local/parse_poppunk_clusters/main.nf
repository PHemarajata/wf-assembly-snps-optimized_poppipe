process PARSE_POPPUNK_CLUSTERS {
    tag "PopPUNK cluster parsing"
    label "process_single"

    conda "conda-forge::python=3.9 conda-forge::pandas=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.3' :
        'staphb/pandas:2.2.3' }"

    input:
    path poppunk_assignments
    val input_file_paths

    output:
    path "cluster_*.txt", emit: cluster_files
    path "cluster_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

import pandas as pd
import os
import glob
from pathlib import Path

# Read PopPUNK assignments
df = pd.read_csv('${poppunk_assignments}')

# Get list of input files from the input_file_paths parameter
input_files = []
file_paths_str = '''${input_file_paths}'''

# Parse the file paths (they are joined with newlines)
for file_path in file_paths_str.strip().split('\\n'):
    if file_path and os.path.isfile(file_path):
        # Check if it's a FASTA file
        for ext in ['.fasta', '.fas', '.fna', '.fsa', '.fa', '.fasta.gz', '.fas.gz', '.fna.gz', '.fsa.gz', '.fa.gz']:
            if file_path.endswith(ext):
                input_files.append(file_path)
                break

# Create mapping of sample names to files with flexible matching
sample_to_file = {}
file_variants = {}

for file_path in input_files:
    file = os.path.basename(file_path)
    # Extract sample name from filename (remove extensions)
    sample_name = Path(file).stem
    
    # Remove common genomic file extensions
    for ext in ['.genomic', '.fasta', '.fas', '.fna', '.fsa', '.fa']:
        if sample_name.endswith(ext):
            sample_name = sample_name[:-len(ext)]
    
    # Create multiple variants of the sample name for matching
    variants = set()
    
    # 1. Original name
    variants.add(sample_name)
    
    # 2. Remove common suffixes
    clean_name = sample_name.replace('-SPAdes', '').replace('_contigs', '').replace('_genomic', '')
    variants.add(clean_name)
    
    # 3. Extract base sample ID (everything before first underscore)
    if '_' in clean_name:
        base_name = clean_name.split('_')[0]
        if len(base_name) > 3:  # Only if meaningful length
            variants.add(base_name)
    
    # 4. Extract base sample ID (everything before first dash)
    if '-' in clean_name:
        base_name = clean_name.split('-')[0]
        if len(base_name) > 3:  # Only if meaningful length
            variants.add(base_name)
    
    # 5. For GCA/GCF accessions, try without version
    if clean_name.startswith(('GCA_', 'GCF_')):
        # Remove version suffix (e.g., GCA_002587385.1_ASM258738v1 -> GCA_002587385)
        parts = clean_name.split('_')
        if len(parts) >= 2:
            accession = '_'.join(parts[:2])  # GCA_002587385
            variants.add(accession)
            # Also try with just the number part
            if '.' in parts[1]:
                number_part = parts[1].split('.')[0]
                variants.add(f"{parts[0]}_{number_part}")
    
    # Store all variants
    for variant in variants:
        if variant:  # Only non-empty variants
            sample_to_file[variant] = file_path
            file_variants[variant] = file_path

# Debug: Print available files and sample mapping
print(f"Available input files: {input_files}")
print(f"Sample to file mapping keys: {list(sample_to_file.keys())}")
print(f"PopPUNK samples: {df['Taxon'].tolist()}")

# Group by cluster
clusters = df.groupby('Cluster')

cluster_summary = []

for cluster_id, group in clusters:
    cluster_files = []
    cluster_filename = f"cluster_{cluster_id}.txt"
    
    with open(cluster_filename, 'w') as f:
        for _, row in group.iterrows():
            sample_name = row['Taxon']  # Adjust column name if different
            matched_file = None
            
            # Try exact match first
            if sample_name in sample_to_file:
                matched_file = sample_to_file[sample_name]
            else:
                # Try fuzzy matching - look for partial matches
                for variant, file_path in sample_to_file.items():
                    if (sample_name in variant or variant in sample_name or 
                        sample_name.replace('_', '').replace('-', '') == variant.replace('_', '').replace('-', '')):
                        matched_file = file_path
                        break
            
            if matched_file:
                f.write(f"{matched_file}\\n")
                cluster_files.append(matched_file)
                print(f"Matched '{sample_name}' to '{matched_file}'")
            else:
                print(f"No match found for sample: '{sample_name}'")
    
    cluster_summary.append({
        'Cluster': cluster_id,
        'Sample_Count': len(cluster_files),
        'Samples': ';'.join([Path(f).stem for f in cluster_files])
    })

# Write summary
summary_df = pd.DataFrame(cluster_summary)
summary_df.to_csv('cluster_summary.tsv', sep='\\t', index=False)

print(f"Created {len(clusters)} cluster files")
for _, row in summary_df.iterrows():
    print(f"Cluster {row['Cluster']}: {row['Sample_Count']} samples")

# Write versions
with open('versions.yml', 'w') as f:
    f.write('"${task.process}":\\n')
    f.write('    python: "3.9"\\n')
    f.write('    pandas: "2.2.3"\\n')
    """
}