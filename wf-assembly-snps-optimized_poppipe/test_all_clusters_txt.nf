#!/usr/bin/env nextflow

/*
 * Test script to verify all_clusters.txt parsing
 * This script tests that the pipeline correctly prioritizes all_clusters.txt over CSV files
 */

nextflow.enable.dsl = 2

// Test parameters
params.outdir = "test_results"

// Create test data with both all_clusters.txt and CSV files
process CREATE_TEST_ALL_CLUSTERS_DATA {
    output:
    path "test_poppipe_output", emit: poppipe_dir
    path "test_assemblies", emit: assembly_dir
    path "rfiles.txt", emit: rfiles

    script:
    """
    # Create test PopPIPE output structure
    mkdir -p test_poppipe_output/strain/cluster_1
    mkdir -p test_poppipe_output/strain/cluster_2
    mkdir -p test_poppipe_output/strain/other
    
    # Create all_clusters.txt (preferred format) - CSV with proper headers
    cat > test_poppipe_output/all_clusters.txt << EOF
Taxon,Strain,Subcluster_1,Subcluster_2,Subcluster_3
sample1,1,1_1,1_1_1,1_1_1_1
sample2,1,1_1,1_1_2,1_1_2_1
sample3,1,1_2,1_2_1,1_2_1_1
sample4,2,2_1,2_1_1,2_1_1_1
sample5,2,2_1,2_1_2,2_1_2_1
sample6,2,2_2,2_2_1,2_2_1_1
EOF

    # Also create a CSV file to test prioritization
    cat > test_poppipe_output/combined_clusters.csv << EOF
Taxon,Cluster
sample1,1
sample2,1
sample3,1
sample4,2
sample5,2
EOF

    # Create test assembly files
    mkdir -p test_assemblies
    
    # Create dummy FASTA files
    for i in {1..6}; do
        cat > test_assemblies/sample\${i}.fasta << EOF
>sample\${i}_contig1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sample\${i}_contig2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
EOF
    done
    
    # Create rfiles.txt
    for i in {1..6}; do
        echo -e "sample\${i}\t\$PWD/test_assemblies/sample\${i}.fasta"
    done > rfiles.txt
    
    # Create test NJ trees in strain directories
    cat > test_poppipe_output/strain/cluster_1/njtree.nwk << 'EOF'
(sample1:0.1,sample2:0.1,(sample3:0.05):0.05);
EOF

    cat > test_poppipe_output/strain/cluster_2/njtree.nwk << 'EOF'
(sample4:0.1,sample5:0.1,sample6:0.1);
EOF

    cat > test_poppipe_output/strain/other/njtree.nwk << 'EOF'
(((sample1:0.1,sample2:0.1):0.05,sample3:0.05):0.1,((sample4:0.1,sample5:0.1,sample6:0.1):0.1):0.1);
EOF
    """
}

// Test the PopPIPE output parser with all_clusters.txt
process TEST_ALL_CLUSTERS_PARSING {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path poppipe_dir
    path rfiles

    output:
    path "cluster_*.txt"
    path "cluster_summary.tsv"
    path "cluster_info.json"
    path "parsing_log.txt"

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
    import json
    import glob
    from pathlib import Path

    # Redirect stdout to log file
    import sys
    log_file = open('parsing_log.txt', 'w')
    sys.stdout = log_file

    print("=== Testing all_clusters.txt parsing ===")
    
    # Find the cluster file with cluster assignments in PopPIPE output
    cluster_files = []
    poppipe_dir = '${poppipe_dir}'

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
            print(f"Found files with pattern '{pattern}': {found_files}")
            break  # Use the first pattern that finds files

    if not cluster_files:
        raise FileNotFoundError(f"No cluster file found in PopPIPE output directory: {poppipe_dir}")

    cluster_file = cluster_files[0]
    print(f"Using cluster file: {cluster_file}")
    
    # Verify that all_clusters.txt was selected over CSV
    if 'all_clusters.txt' in cluster_file:
        print("SUCCESS: all_clusters.txt was correctly prioritized!")
    else:
        print(f"WARNING: Expected all_clusters.txt but got: {cluster_file}")

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
                df = pd.read_csv(cluster_file, sep='\\t')
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
    print(f"Data preview:")
    print(df.head())

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
    cluster_candidates = ['Cluster', 'cluster', 'Cluster_ID', 'cluster_id', 'strain', 'Strain', 'group', 'Group']
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

    # Read rfiles.txt for sample-to-file mapping
    rfiles_path = '${rfiles}'
    sample_to_file = {}
    
    print(f"Reading sample-to-file mapping from: {rfiles_path}")
    with open(rfiles_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\\t')
            if len(parts) >= 2:
                sample_name = parts[0].strip()
                file_path = parts[1].strip()
                
                if os.path.isfile(file_path):
                    sample_to_file[sample_name] = file_path
                    print(f"Mapped '{sample_name}' to '{file_path}'")
                else:
                    print(f"Warning: File not found for sample '{sample_name}': {file_path}")

    print(f"Loaded {len(sample_to_file)} sample-to-file mappings")

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
                
                # Try exact match
                if sample_name in sample_to_file:
                    matched_file = sample_to_file[sample_name]
                
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

    # Write cluster info JSON
    cluster_info = {"test": "all_clusters.txt parsing"}
    with open('cluster_info.json', 'w') as f:
        json.dump(cluster_info, f, indent=2)

    print(f"\\nCreated {len(clusters)} cluster files:")
    for _, row in summary_df.iterrows():
        print(f"Cluster {row['Cluster']}: {row['Sample_Count']} samples")
    
    # Verify expected results
    expected_clusters = {1: 3, 2: 3}  # Based on all_clusters.txt (not CSV)
    actual_clusters = {int(row['Cluster']): row['Sample_Count'] for _, row in summary_df.iterrows()}
    
    print(f"\\nExpected clusters: {expected_clusters}")
    print(f"Actual clusters: {actual_clusters}")
    
    if actual_clusters == expected_clusters:
        print("SUCCESS: Cluster parsing matches expected results from all_clusters.txt!")
    else:
        print("WARNING: Cluster results don't match expected values")
    
    print("\\n=== Test completed ===")
    
    log_file.close()
    """
}

workflow TEST_ALL_CLUSTERS_TXT {
    
    // Create test data
    CREATE_TEST_ALL_CLUSTERS_DATA()
    
    // Test the parser
    TEST_ALL_CLUSTERS_PARSING(
        CREATE_TEST_ALL_CLUSTERS_DATA.out.poppipe_dir,
        CREATE_TEST_ALL_CLUSTERS_DATA.out.rfiles
    )
    
    // Display results
    TEST_ALL_CLUSTERS_PARSING.out[0].view { "Cluster files: $it" }
    TEST_ALL_CLUSTERS_PARSING.out[3].view { "Log content:\n" + it.text }
}

workflow {
    TEST_ALL_CLUSTERS_TXT()
}