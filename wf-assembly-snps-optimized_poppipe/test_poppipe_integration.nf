#!/usr/bin/env nextflow

/*
 * Test script for PopPIPE-bp integration
 * This script validates the new modules and workflow
 */

nextflow.enable.dsl = 2

// Test parameters
params.test_poppipe_output = "test_poppipe_output"
params.test_assemblies = "test_assemblies"
params.outdir = "test_results"

// Create test data
process CREATE_TEST_DATA {
    output:
    path "test_poppipe_output", emit: poppipe_dir
    path "test_assemblies", emit: assembly_dir
    path "assembly_list.txt", emit: assembly_list

    script:
    """
    # Create test PopPIPE output structure
    mkdir -p test_poppipe_output/strain/cluster_1
    mkdir -p test_poppipe_output/strain/cluster_2
    mkdir -p test_poppipe_output/strain/other
    
    # Create test CSV file
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
    
    # Create assembly list
    ls -1 \$PWD/test_assemblies/*.fasta > assembly_list.txt
    
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
    """
}

// Test the PopPIPE output parser
process TEST_PARSE_POPPIPE_OUTPUT {
    input:
    path poppipe_dir
    val assembly_paths

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

    # Find CSV file
    csv_files = glob.glob('${poppipe_dir}/*.csv')
    if not csv_files:
        raise FileNotFoundError("No CSV file found")
    
    csv_file = csv_files[0]
    print(f"Using CSV file: {csv_file}")
    
    # Read cluster assignments
    df = pd.read_csv(csv_file)
    print(f"CSV shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")
    
    # Parse assembly paths
    assembly_files = '${assembly_paths}'.strip().split()
    print(f"Assembly files: {assembly_files}")
    
    # Create sample to file mapping
    sample_to_file = {}
    for file_path in assembly_files:
        if os.path.isfile(file_path):
            sample_name = Path(file_path).stem
            sample_to_file[sample_name] = file_path
    
    print(f"Sample mapping: {sample_to_file}")
    
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
        
        cluster_summary.append({
            'Cluster': cluster_id,
            'Sample_Count': len(cluster_files),
            'Samples': ';'.join([Path(f).stem for f in cluster_files])
        })
    
    # Write summary
    summary_df = pd.DataFrame(cluster_summary)
    summary_df.to_csv('cluster_summary.tsv', sep='\\t', index=False)
    
    # Write cluster info
    cluster_info = {"test": "data"}
    with open('cluster_info.json', 'w') as f:
        json.dump(cluster_info, f, indent=2)
    
    print("Test parsing completed successfully")
    """
}

workflow TEST_POPPIPE_INTEGRATION {
    
    // Create test data
    CREATE_TEST_DATA()
    
    // Read assembly list and convert to space-separated string
    assembly_paths = CREATE_TEST_DATA.out.assembly_list
        .splitText()
        .map { it.trim() }
        .collect()
        .map { it.join(' ') }
    
    // Test the parser
    TEST_PARSE_POPPIPE_OUTPUT(
        CREATE_TEST_DATA.out.poppipe_dir,
        assembly_paths
    )
    
    // Display results
    TEST_PARSE_POPPIPE_OUTPUT.out[0].view { "Cluster files: $it" }
    TEST_PARSE_POPPIPE_OUTPUT.out[1].view { "Summary: " + it.text }
}

workflow {
    TEST_POPPIPE_INTEGRATION()
}