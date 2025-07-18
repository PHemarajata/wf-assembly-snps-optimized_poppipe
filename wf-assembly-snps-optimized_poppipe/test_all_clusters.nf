#!/usr/bin/env nextflow

/*
 * Test script for all_clusters.txt format validation
 * Tests the PopPIPE-bp integration with all_clusters.txt format
 */

nextflow.enable.dsl = 2

// Test parameters
params.outdir = "test_results_all_clusters"

// Create test data with all_clusters.txt format
process CREATE_ALL_CLUSTERS_TEST_DATA {
    output:
    path "test_poppipe_output", emit: poppipe_dir
    path "test_assemblies", emit: assembly_dir
    path "rfiles.txt", emit: rfiles

    script:
    """
    # Create test PopPIPE output structure
    mkdir -p test_poppipe_output/strains/1
    mkdir -p test_poppipe_output/strains/2
    mkdir -p test_poppipe_output/strains/other
    
    # Create all_clusters.txt with the specified format
    cat > test_poppipe_output/all_clusters.txt << EOF
Taxon,Strain,Subcluster_1,Subcluster_2,Subcluster_3
sample001,1,1_1,1_1_1,1_1_1_1
sample002,1,1_1,1_1_1,1_1_1_2
sample003,1,1_2,1_2_1,1_2_1_1
sample004,2,2_1,2_1_1,2_1_1_1
sample005,2,2_1,2_1_2,2_1_2_1
sample006,2,2_2,2_2_1,2_2_1_1
EOF

    # Create test assembly files
    mkdir -p test_assemblies
    
    # Create dummy FASTA files
    for i in {1..6}; do
        sample_name="sample00\${i}"
        cat > test_assemblies/\${sample_name}.fasta << EOF
>\${sample_name}_contig1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>\${sample_name}_contig2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
EOF
    done
    
    # Create rfiles.txt for TSV input testing
    cat > rfiles.txt << EOF
sample001	\$PWD/test_assemblies/sample001.fasta
sample002	\$PWD/test_assemblies/sample002.fasta
sample003	\$PWD/test_assemblies/sample003.fasta
sample004	\$PWD/test_assemblies/sample004.fasta
sample005	\$PWD/test_assemblies/sample005.fasta
sample006	\$PWD/test_assemblies/sample006.fasta
EOF
    
    # Create PopPIPE strain directory structure with expected files
    # Strain 1 files
    cat > test_poppipe_output/strains/1/njtree.nwk << 'EOF'
(sample001:0.1,sample002:0.1,(sample003:0.05):0.05);
EOF

    cat > test_poppipe_output/strains/1/besttree.nwk << 'EOF'
(sample001:0.12,sample002:0.11,(sample003:0.06):0.04);
EOF

    cat > test_poppipe_output/strains/1/rfile.txt << 'EOF'
sample001	/path/to/sample001.fasta
sample002	/path/to/sample002.fasta
sample003	/path/to/sample003.fasta
EOF

    cat > test_poppipe_output/strains/1/names.txt << 'EOF'
sample001
sample002
sample003
EOF

    cat > test_poppipe_output/strains/1/fastbaps_clusters.txt << 'EOF'
sample001	1
sample002	1
sample003	2
EOF

    # Create dummy files for strain 1
    touch test_poppipe_output/strains/1/align_variants.aln
    touch test_poppipe_output/strains/1/dists.npy
    touch test_poppipe_output/strains/1/dists.pkl
    touch test_poppipe_output/strains/1/split_kmers.skf

    # Strain 2 files
    cat > test_poppipe_output/strains/2/njtree.nwk << 'EOF'
(sample004:0.1,(sample005:0.08,sample006:0.09):0.05);
EOF

    cat > test_poppipe_output/strains/2/besttree.nwk << 'EOF'
(sample004:0.11,(sample005:0.09,sample006:0.10):0.06);
EOF

    cat > test_poppipe_output/strains/2/rfile.txt << 'EOF'
sample004	/path/to/sample004.fasta
sample005	/path/to/sample005.fasta
sample006	/path/to/sample006.fasta
EOF

    cat > test_poppipe_output/strains/2/names.txt << 'EOF'
sample004
sample005
sample006
EOF

    cat > test_poppipe_output/strains/2/fastbaps_clusters.txt << 'EOF'
sample004	1
sample005	2
sample006	2
EOF

    # Create dummy files for strain 2
    touch test_poppipe_output/strains/2/align_variants.aln
    touch test_poppipe_output/strains/2/dists.npy
    touch test_poppipe_output/strains/2/dists.pkl
    touch test_poppipe_output/strains/2/split_kmers.skf

    # Create backbone tree
    cat > test_poppipe_output/strains/other/njtree.nwk << 'EOF'
(((sample001:0.1,sample002:0.1):0.05,sample003:0.05):0.1,((sample004:0.1,(sample005:0.08,sample006:0.09):0.05):0.1):0.1);
EOF

    echo "Created test data with all_clusters.txt format"
    """
}

// Test the PopPIPE output parser with all_clusters.txt
process TEST_ALL_CLUSTERS_PARSING {
    input:
    path poppipe_dir
    path rfiles

    output:
    path "cluster_*.txt"
    path "cluster_summary.tsv"
    path "cluster_info.json"

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
    import json
    import glob
    from pathlib import Path

    # Find the cluster file - prioritize all_clusters.txt
    poppipe_dir = '${poppipe_dir}'
    cluster_files = []
    
    # Look for all_clusters.txt first
    search_patterns = [
        os.path.join(poppipe_dir, 'all_clusters.txt'),
        os.path.join(poppipe_dir, 'output', 'all_clusters.txt'),
        os.path.join(poppipe_dir, '*_clusters.csv'),
        os.path.join(poppipe_dir, '*.csv')
    ]

    for pattern in search_patterns:
        found_files = glob.glob(pattern)
        if found_files:
            cluster_files.extend(found_files)
            break

    if not cluster_files:
        raise FileNotFoundError(f"No cluster file found in PopPIPE output directory: {poppipe_dir}")

    cluster_file = cluster_files[0]
    print(f"Using cluster file: {cluster_file}")

    # Read the all_clusters.txt file
    df = pd.read_csv(cluster_file)
    print(f"Cluster file columns: {df.columns.tolist()}")
    print(f"Cluster file shape: {df.shape}")
    print(f"First few rows:")
    print(df.head())

    # Identify columns - for all_clusters.txt format
    sample_col = 'Taxon'  # Standard column name
    cluster_col = 'Strain'  # Primary cluster column in all_clusters.txt
    
    if sample_col not in df.columns:
        raise ValueError(f"Expected 'Taxon' column not found. Available columns: {df.columns.tolist()}")
    if cluster_col not in df.columns:
        raise ValueError(f"Expected 'Strain' column not found. Available columns: {df.columns.tolist()}")

    print(f"Using sample column: {sample_col}, cluster column: {cluster_col}")

    # Read TSV file for sample-to-file mapping
    rfiles_path = '${rfiles}'
    sample_to_file = {}
    
    with open(rfiles_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) >= 2:
                sample_name = parts[0].strip()
                file_path = parts[1].strip()
                
                if os.path.isfile(file_path):
                    sample_to_file[sample_name] = file_path
                    print(f"Mapped '{sample_name}' to '{file_path}'")
                else:
                    print(f"Warning: File not found for sample '{sample_name}': {file_path}")

    print(f"Loaded {len(sample_to_file)} sample-to-file mappings from TSV")

    # Find strain directories
    strain_paths = [
        os.path.join(poppipe_dir, 'strain'),
        os.path.join(poppipe_dir, 'strains'),
        os.path.join(poppipe_dir, 'output', 'strain'),
        os.path.join(poppipe_dir, 'output', 'strains')
    ]

    cluster_info = {}
    for strain_path in strain_paths:
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

    # Group by strain (cluster) and create cluster files
    clusters = df.groupby(cluster_col)
    cluster_summary = []

    for cluster_id, group in clusters:
        cluster_files = []
        cluster_filename = f"cluster_{cluster_id}.txt"
        
        print(f"\\nProcessing cluster {cluster_id} with {len(group)} samples:")
        
        with open(cluster_filename, 'w') as f:
            for _, row in group.iterrows():
                sample_name = row[sample_col]
                print(f"  Sample: {sample_name}")
                
                # Show subcluster information if available
                if 'Subcluster_1' in df.columns:
                    print(f"    Subcluster_1: {row.get('Subcluster_1', 'N/A')}")
                if 'Subcluster_2' in df.columns:
                    print(f"    Subcluster_2: {row.get('Subcluster_2', 'N/A')}")
                if 'Subcluster_3' in df.columns:
                    print(f"    Subcluster_3: {row.get('Subcluster_3', 'N/A')}")
                
                matched_file = sample_to_file.get(sample_name)
                if matched_file:
                    f.write(f"{matched_file}\\n")
                    cluster_files.append(matched_file)
                    print(f"    Matched to: {matched_file}")
                else:
                    print(f"    No match found for sample: {sample_name}")
        
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

    print(f"\\nCreated {len(clusters)} cluster files")
    for _, row in summary_df.iterrows():
        print(f"Cluster {row['Cluster']}: {row['Sample_Count']} samples, PopPIPE data: {row['Has_PopPIPE_Data']}")

    print("\\nall_clusters.txt parsing test completed successfully")
    """
}

workflow TEST_ALL_CLUSTERS {
    
    // Create test data
    CREATE_ALL_CLUSTERS_TEST_DATA()
    
    // Test the parser
    TEST_ALL_CLUSTERS_PARSING(
        CREATE_ALL_CLUSTERS_TEST_DATA.out.poppipe_dir,
        CREATE_ALL_CLUSTERS_TEST_DATA.out.rfiles
    )
    
    // Display results
    TEST_ALL_CLUSTERS_PARSING.out[0].view { "Cluster files: $it" }
    TEST_ALL_CLUSTERS_PARSING.out[1].view { "Summary file: $it" }
}

workflow {
    TEST_ALL_CLUSTERS()
}