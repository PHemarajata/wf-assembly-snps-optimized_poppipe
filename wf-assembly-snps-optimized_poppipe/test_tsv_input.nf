#!/usr/bin/env nextflow

/*
 * Test script for PopPIPE-bp TSV input functionality
 * This script validates the TSV file input capability
 */

nextflow.enable.dsl = 2

// Test parameters
params.outdir = "test_tsv_results"

// Create test data with TSV input
process CREATE_TEST_TSV_DATA {
    output:
    path "test_poppipe_output", emit: poppipe_dir
    path "test_assemblies", emit: assembly_dir
    path "rfiles.txt", emit: tsv_file

    script:
    """
    # Create test PopPIPE output structure
    mkdir -p test_poppipe_output/strain/cluster_1
    mkdir -p test_poppipe_output/strain/cluster_2
    mkdir -p test_poppipe_output/strain/other
    
    # Create test CSV file with cluster assignments
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
    for i in {1..5}; do
        cat > test_assemblies/sample\${i}.fasta << EOF
>sample\${i}_contig1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sample\${i}_contig2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
EOF
    done
    
    # Create TSV file (rfiles.txt format) with sample-to-file mapping
    cat > rfiles.txt << EOF
sample1	\$PWD/test_assemblies/sample1.fasta
sample2	\$PWD/test_assemblies/sample2.fasta
sample3	\$PWD/test_assemblies/sample3.fasta
sample4	\$PWD/test_assemblies/sample4.fasta
sample5	\$PWD/test_assemblies/sample5.fasta
EOF
    
    # Create test NJ trees in strain directories
    cat > test_poppipe_output/strain/cluster_1/njtree.nwk << 'EOF'
(sample1:0.1,sample2:0.1,(sample3:0.05):0.05);
EOF

    cat > test_poppipe_output/strain/cluster_2/njtree.nwk << 'EOF'
(sample4:0.1,sample5:0.1);
EOF

    cat > test_poppipe_output/strain/other/njtree.nwk << 'EOF'
(((sample1:0.1,sample2:0.1):0.05,sample3:0.05):0.1,((sample4:0.1,sample5:0.1):0.1):0.1);
EOF

    echo "Created test data with TSV input file"
    """
}

// Test the PopPIPE output parser with TSV input
process TEST_TSV_PARSING {
    input:
    path poppipe_dir
    path tsv_file

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

    # Log file for debugging
    log_file = open('parsing_log.txt', 'w')
    
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
        
        # Read TSV file
        tsv_file_path = '${tsv_file}'
        log_print(f"Reading TSV file: {tsv_file_path}")
        
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
                    
                    if os.path.isfile(file_path):
                        sample_to_file[sample_name] = file_path
                        log_print(f"Mapped '{sample_name}' to '{file_path}'")
                    else:
                        log_print(f"Warning: File not found for sample '{sample_name}': {file_path}")
                else:
                    log_print(f"Warning: Invalid format at line {line_num}: {line}")
        
        log_print(f"Loaded {len(sample_to_file)} sample-to-file mappings from TSV")
        
        # Group by cluster
        clusters = df.groupby('Cluster')
        cluster_summary = []
        
        for cluster_id, group in clusters:
            cluster_files = []
            cluster_filename = f"cluster_{cluster_id}.txt"
            
            with open(cluster_filename, 'w') as f:
                for _, row in group.iterrows():
                    sample_name = row['Taxon']
                    if sample_name in sample_to_file:
                        matched_file = sample_to_file[sample_name]
                        f.write(f"{matched_file}\\n")
                        cluster_files.append(matched_file)
                        log_print(f"Matched '{sample_name}' to '{matched_file}'")
                    else:
                        log_print(f"No match found for sample: '{sample_name}'")
            
            cluster_summary.append({
                'Cluster': cluster_id,
                'Sample_Count': len(cluster_files),
                'Samples': ';'.join([Path(f).stem for f in cluster_files])
            })
        
        # Write summary
        summary_df = pd.DataFrame(cluster_summary)
        summary_df.to_csv('cluster_summary.tsv', sep='\\t', index=False)
        
        # Write cluster info
        cluster_info = {"test": "TSV input successful"}
        with open('cluster_info.json', 'w') as f:
            json.dump(cluster_info, f, indent=2)
        
        log_print("TSV parsing test completed successfully")
        
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

workflow TEST_TSV_INPUT {
    
    // Create test data
    CREATE_TEST_TSV_DATA()
    
    // Test TSV parsing
    TEST_TSV_PARSING(
        CREATE_TEST_TSV_DATA.out.poppipe_dir,
        CREATE_TEST_TSV_DATA.out.tsv_file
    )
    
    // Display results
    TEST_TSV_PARSING.out[0].view { "Cluster files: $it" }
    TEST_TSV_PARSING.out[1].view { "Summary file: $it" }
    TEST_TSV_PARSING.out[3].view { "Log: " + it.text }
}

workflow {
    TEST_TSV_INPUT()
}