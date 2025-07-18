process PARSE_POPPIPE_OUTPUT {
    tag "PopPIPE output parsing"
    label "process_single"

    conda "conda-forge::python=3.9 conda-forge::pandas=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.3' :
        'staphb/pandas:2.2.3' }"

    input:
    path poppipe_output_dir
    path input_source  // Can be either a TSV file (rfiles.txt) or directory

    output:
    path "cluster_*.txt", emit: cluster_files
    path "cluster_summary.tsv", emit: summary
    path "cluster_info.json", emit: cluster_info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

import pandas as pd
import os
import json
import glob
from pathlib import Path

# Find the cluster file with cluster assignments in PopPIPE output
cluster_files = []
poppipe_dir = '${poppipe_output_dir}'

# Look for cluster files in priority order - all_clusters.txt first, then CSV files
search_patterns = [
    # Priority 1: all_clusters.txt (PopPIPE standard output)
    os.path.join(poppipe_dir, 'all_clusters.txt'),
    os.path.join(poppipe_dir, 'output', 'all_clusters.txt'),
    # Priority 2: Other cluster files
    os.path.join(poppipe_dir, '*_clusters.csv'),
    os.path.join(poppipe_dir, '*.csv'),
    os.path.join(poppipe_dir, 'output', '*_clusters.csv'),
    os.path.join(poppipe_dir, 'output', '*.csv'),
    # Priority 3: TXT files
    os.path.join(poppipe_dir, '*_clusters.txt'),
    os.path.join(poppipe_dir, 'output', '*_clusters.txt')
]

for pattern in search_patterns:
    found_files = glob.glob(pattern)
    if found_files:
        cluster_files.extend(found_files)
        break  # Use the first pattern that finds files

if not cluster_files:
    raise FileNotFoundError(f"No cluster file found in PopPIPE output directory: {poppipe_dir}")

cluster_file = cluster_files[0]
print(f"Using cluster file: {cluster_file}")

# Determine file format and read cluster assignments
file_ext = Path(cluster_file).suffix.lower()
if file_ext == '.txt':
    # Handle TXT format - could be tab-separated or other format
    try:
        # Try reading as CSV first (comma-separated)
        df = pd.read_csv(cluster_file)
        print(f"Successfully read TXT file as CSV format")
    except:
        try:
            # Try reading as TSV (tab-separated)
            df = pd.read_csv(cluster_file, sep='\t')
            print(f"Successfully read TXT file as TSV format")
        except:
            # Try reading with various separators
            for sep in [' ', '|', ';']:
                try:
                    df = pd.read_csv(cluster_file, sep=sep)
                    print(f"Successfully read TXT file with separator '{sep}'")
                    break
                except:
                    continue
            else:
                raise ValueError(f"Could not parse TXT file: {cluster_file}")
else:
    # Handle CSV format
    df = pd.read_csv(cluster_file)
print(f"Cluster file columns: {df.columns.tolist()}")
print(f"Cluster file shape: {df.shape}")

# Try to identify the correct columns
sample_col = None
cluster_col = None

# Common column names for samples (expanded for all_clusters.txt format)
sample_candidates = ['Taxon', 'Sample', 'sample', 'ID', 'id', 'Name', 'name', 'sample_id', 'Sample_ID']
for col in sample_candidates:
    if col in df.columns:
        sample_col = col
        break

# Common column names for clusters (expanded for all_clusters.txt format)
# Note: 'Strain' is the primary cluster column in all_clusters.txt
cluster_candidates = ['Strain', 'strain', 'Cluster', 'cluster', 'Cluster_ID', 'cluster_id', 'group', 'Group']
for col in cluster_candidates:
    if col in df.columns:
        cluster_col = col
        break

# If standard column names not found, try positional detection
if not sample_col and len(df.columns) >= 2:
    # Assume first column is sample, second is cluster
    sample_col = df.columns[0]
    print(f"Using first column as sample column: {sample_col}")

if not cluster_col and len(df.columns) >= 2:
    # Assume second column is cluster
    cluster_col = df.columns[1]
    print(f"Using second column as cluster column: {cluster_col}")

if not sample_col or not cluster_col:
    print(f"Available columns: {df.columns.tolist()}")
    raise ValueError(f"Could not identify sample column ({sample_col}) or cluster column ({cluster_col})")

print(f"Using sample column: {sample_col}, cluster column: {cluster_col}")

# Determine if input_source is a TSV file or directory
input_source_path = '${input_source}'
sample_to_file = {}

if os.path.isfile(input_source_path):
    # Input is a TSV file (rfiles.txt, combined.txt, etc.)
    print(f"Reading sample-to-file mapping from TSV: {input_source_path}")
    
    # Get the directory containing the TSV file for resolving relative paths
    tsv_dir = os.path.dirname(os.path.abspath(input_source_path))
    print(f"TSV file directory: {tsv_dir}")
    
    try:
        # Read TSV file - assume tab-separated with sample_name\tfile_path format
        with open(input_source_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    sample_name = parts[0].strip()
                    file_path = parts[1].strip()
                    
                    # If file_path is relative, make it relative to the TSV file directory
                    if not os.path.isabs(file_path):
                        file_path = os.path.join(tsv_dir, file_path)
                    
                    # Normalize the path to resolve any .. or . components
                    file_path = os.path.normpath(file_path)
                    
                    # Verify file exists
                    if os.path.isfile(file_path):
                        sample_to_file[sample_name] = file_path
                        print(f"Mapped '{sample_name}' to '{file_path}'")
                    else:
                        print(f"Warning: File not found for sample '{sample_name}': {file_path}")
                else:
                    print(f"Warning: Invalid format at line {line_num}: {line}")
        
        print(f"Loaded {len(sample_to_file)} sample-to-file mappings from TSV")
        
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        raise

elif os.path.isdir(input_source_path):
    # Input is a directory - use original file discovery logic
    print(f"Discovering assembly files in directory: {input_source_path}")
    
    input_files = []
    for ext in ['.fasta', '.fas', '.fna', '.fsa', '.fa', '.fasta.gz', '.fas.gz', '.fna.gz', '.fsa.gz', '.fa.gz']:
        pattern = os.path.join(input_source_path, f"*{ext}")
        input_files.extend(glob.glob(pattern))
    
    print(f"Found {len(input_files)} input assembly files")
    
    # Create mapping of sample names to files using filename-based matching
    for file_path in input_files:
        file = os.path.basename(file_path)
        sample_name = Path(file).stem
        
        # Remove common genomic file extensions
        for ext in ['.genomic', '.fasta', '.fas', '.fna', '.fsa', '.fa']:
            if sample_name.endswith(ext):
                sample_name = sample_name[:-len(ext)]
        
        # Create multiple variants for matching
        variants = set([sample_name])
        
        # Remove common suffixes
        clean_name = sample_name.replace('-SPAdes', '').replace('_contigs', '').replace('_genomic', '')
        variants.add(clean_name)
        
        # Extract base sample ID
        if '_' in clean_name:
            base_name = clean_name.split('_')[0]
            if len(base_name) > 3:
                variants.add(base_name)
        
        if '-' in clean_name:
            base_name = clean_name.split('-')[0]
            if len(base_name) > 3:
                variants.add(base_name)
        
        # For GCA/GCF accessions
        if clean_name.startswith(('GCA_', 'GCF_')):
            parts = clean_name.split('_')
            if len(parts) >= 2:
                accession = '_'.join(parts[:2])
                variants.add(accession)
        
        for variant in variants:
            if variant:
                sample_to_file[variant] = file_path

else:
    raise ValueError(f"Input source is neither a file nor directory: {input_source_path}")

# Find strain directories in PopPIPE output
strain_dirs = []
possible_strain_paths = [
    os.path.join(poppipe_dir, 'strain'),
    os.path.join(poppipe_dir, 'strains'),
    os.path.join(poppipe_dir, 'output', 'strain'),
    os.path.join(poppipe_dir, 'output', 'strains')
]

cluster_info = {}
for strain_path in possible_strain_paths:
    if os.path.exists(strain_path):
        for cluster_dir in os.listdir(strain_path):
            cluster_path = os.path.join(strain_path, cluster_dir)
            if os.path.isdir(cluster_path):
                cluster_info[cluster_dir] = {
                    'path': cluster_path,
                    'files': os.listdir(cluster_path)
                }
        break

print(f"Found {len(cluster_info)} cluster directories in strain folder")

# Group by cluster and create cluster files
clusters = df.groupby(cluster_col)
cluster_summary = []

for cluster_id, group in clusters:
    cluster_files = []
    cluster_filename = f"cluster_{cluster_id}.txt"
    
    with open(cluster_filename, 'w') as f:
        for _, row in group.iterrows():
            sample_name = row[sample_col]
            matched_file = None
            
            # Try exact match first
            if sample_name in sample_to_file:
                matched_file = sample_to_file[sample_name]
            else:
                # Try fuzzy matching (only if input was a directory, not TSV)
                if os.path.isdir(input_source_path):
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
        'Samples': ';'.join([Path(f).stem for f in cluster_files]),
        'Has_PopPIPE_Data': str(cluster_id) in cluster_info
    })

# Write summary
summary_df = pd.DataFrame(cluster_summary)
summary_df.to_csv('cluster_summary.tsv', sep='\\t', index=False)

# Write cluster info JSON
with open('cluster_info.json', 'w') as f:
    json.dump(cluster_info, f, indent=2)

print(f"Created {len(clusters)} cluster files")
for _, row in summary_df.iterrows():
    print(f"Cluster {row['Cluster']}: {row['Sample_Count']} samples, PopPIPE data: {row['Has_PopPIPE_Data']}")

# Write versions
with open('versions.yml', 'w') as f:
    f.write('"${task.process}":\\n')
    f.write('    python: "3.9"\\n')
    f.write('    pandas: "2.2.3"\\n')
    """
}