#!/usr/bin/env nextflow

/*
 * Test script for relative path functionality in TSV files
 * This demonstrates how file paths in rfiles.txt/combined_rfile.txt 
 * are resolved relative to the TSV file location
 */

nextflow.enable.dsl = 2

// Test parameters
params.outdir = "test_relative_paths_results"

// Create test data with relative paths in TSV
process CREATE_RELATIVE_PATH_TEST_DATA {
    output:
    path "test_project", emit: project_dir
    path "test_project/data/rfiles.txt", emit: tsv_file
    path "test_project/poppipe_output", emit: poppipe_dir

    script:
    """
    # Create a test project structure with subdirectories
    mkdir -p test_project/data
    mkdir -p test_project/assemblies/batch1
    mkdir -p test_project/assemblies/batch2
    mkdir -p test_project/poppipe_output/strain/cluster_1
    mkdir -p test_project/poppipe_output/strain/cluster_2
    
    # Create test assembly files in different subdirectories
    cat > test_project/assemblies/batch1/sample1.fasta << EOF
>sample1_contig1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

    cat > test_project/assemblies/batch1/sample2.fasta << EOF
>sample2_contig1
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
EOF

    cat > test_project/assemblies/batch2/sample3.fasta << EOF
>sample3_contig1
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
EOF

    cat > test_project/assemblies/batch2/sample4.fasta << EOF
>sample4_contig1
CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG
CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG
EOF

    # Create cluster assignments CSV
    cat > test_project/poppipe_output/combined_clusters.csv << EOF
Taxon,Cluster
sample1,1
sample2,1
sample3,2
sample4,2
EOF

    # Create rfiles.txt with RELATIVE paths (relative to the rfiles.txt location)
    # The rfiles.txt is in test_project/data/
    # The assemblies are in test_project/assemblies/
    # So relative paths should be ../assemblies/...
    cat > test_project/data/rfiles.txt << EOF
sample1	../assemblies/batch1/sample1.fasta
sample2	../assemblies/batch1/sample2.fasta
sample3	../assemblies/batch2/sample3.fasta
sample4	../assemblies/batch2/sample4.fasta
EOF

    echo "Created test project with relative paths in rfiles.txt"
    echo "TSV file location: test_project/data/rfiles.txt"
    echo "Assembly files referenced with relative paths: ../assemblies/..."
    """
}

// Test the relative path resolution
process TEST_RELATIVE_PATH_RESOLUTION {
    input:
    path project_dir
    path tsv_file
    path poppipe_dir

    output:
    path "cluster_*.txt"
    path "cluster_summary.tsv"
    path "relative_path_test_log.txt"

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
    import glob
    from pathlib import Path

    # Log file for debugging
    log_file = open('relative_path_test_log.txt', 'w')
    
    def log_print(msg):
        print(msg)
        log_file.write(msg + '\\n')
        log_file.flush()

    try:
        # Find CSV file
        csv_files = glob.glob('${poppipe_dir}/*.csv')
        if not csv_files:
            raise FileNotFoundError("No CSV file found")
        
        csv_file = csv_files[0]
        log_print(f"Using CSV file: {csv_file}")
        
        # Read cluster assignments
        df = pd.read_csv(csv_file)
        log_print(f"CSV shape: {df.shape}")
        log_print(f"Columns: {df.columns.tolist()}")
        
        # Read TSV file with relative path resolution
        tsv_file_path = '${tsv_file}'
        log_print(f"Reading TSV file: {tsv_file_path}")
        
        # Get the directory containing the TSV file for resolving relative paths
        tsv_dir = os.path.dirname(os.path.abspath(tsv_file_path))
        log_print(f"TSV file directory: {tsv_dir}")
        
        sample_to_file = {}
        with open(tsv_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\\t')
                if len(parts) >= 2:
                    sample_name = parts[0].strip()
                    file_path = parts[1].strip()
                    
                    log_print(f"Original path from TSV: {file_path}")
                    
                    # If file_path is relative, make it relative to the TSV file directory
                    if not os.path.isabs(file_path):
                        resolved_path = os.path.join(tsv_dir, file_path)
                        log_print(f"Resolved relative path: {resolved_path}")
                    else:
                        resolved_path = file_path
                        log_print(f"Absolute path used as-is: {resolved_path}")
                    
                    # Normalize the path to resolve any .. or . components
                    normalized_path = os.path.normpath(resolved_path)
                    log_print(f"Normalized path: {normalized_path}")
                    
                    if os.path.isfile(normalized_path):
                        sample_to_file[sample_name] = normalized_path
                        log_print(f"✓ Successfully mapped '{sample_name}' to '{normalized_path}'")
                    else:
                        log_print(f"✗ File not found for sample '{sample_name}': {normalized_path}")
                else:
                    log_print(f"Warning: Invalid format at line {line_num}: {line}")
        
        log_print(f"\\nLoaded {len(sample_to_file)} sample-to-file mappings from TSV")
        
        # Group by cluster
        clusters = df.groupby('Cluster')
        cluster_summary = []
        
        for cluster_id, group in clusters:
            cluster_files = []
            cluster_filename = f"cluster_{cluster_id}.txt"
            
            log_print(f"\\nProcessing cluster {cluster_id}:")
            with open(cluster_filename, 'w') as f:
                for _, row in group.iterrows():
                    sample_name = row['Taxon']
                    if sample_name in sample_to_file:
                        matched_file = sample_to_file[sample_name]
                        f.write(f"{matched_file}\\n")
                        cluster_files.append(matched_file)
                        log_print(f"  ✓ {sample_name} -> {matched_file}")
                    else:
                        log_print(f"  ✗ No match found for sample: {sample_name}")
            
            cluster_summary.append({
                'Cluster': cluster_id,
                'Sample_Count': len(cluster_files),
                'Samples': ';'.join([Path(f).stem for f in cluster_files])
            })
        
        # Write summary
        summary_df = pd.DataFrame(cluster_summary)
        summary_df.to_csv('cluster_summary.tsv', sep='\\t', index=False)
        
        log_print("\\n=== RELATIVE PATH TEST SUMMARY ===")
        log_print("✓ Relative path resolution test completed successfully")
        
        # Print summary
        for _, row in summary_df.iterrows():
            log_print(f"Cluster {row['Cluster']}: {row['Sample_Count']} samples")
    
    except Exception as e:
        log_print(f"Error: {e}")
        raise
    finally:
        log_file.close()
    """
}

workflow TEST_RELATIVE_PATHS {
    
    // Create test data with relative paths
    CREATE_RELATIVE_PATH_TEST_DATA()
    
    // Test relative path resolution
    TEST_RELATIVE_PATH_RESOLUTION(
        CREATE_RELATIVE_PATH_TEST_DATA.out.project_dir,
        CREATE_RELATIVE_PATH_TEST_DATA.out.tsv_file,
        CREATE_RELATIVE_PATH_TEST_DATA.out.poppipe_dir
    )
    
    // Display results
    TEST_RELATIVE_PATH_RESOLUTION.out[2].view { "Test log:\n" + it.text }
}

workflow {
    TEST_RELATIVE_PATHS()
}