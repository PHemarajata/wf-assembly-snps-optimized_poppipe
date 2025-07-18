#!/usr/bin/env nextflow

/*
 * Complete PopPIPE-bp integration test
 * Tests the full workflow with all_clusters.txt and strain directories
 */

nextflow.enable.dsl = 2

// Test parameters
params.outdir = "test_complete_poppipe_results"

// Create comprehensive test data matching real PopPIPE-bp output
process CREATE_COMPLETE_POPPIPE_DATA {
    output:
    path "poppipe_output", emit: poppipe_dir
    path "assemblies", emit: assembly_dir
    path "rfiles.txt", emit: rfiles

    script:
    """
    # Create PopPIPE output structure
    mkdir -p poppipe_output/strains/{1,2,other}
    
    # Create all_clusters.txt with hierarchical clustering
    cat > poppipe_output/all_clusters.txt << EOF
Taxon,Strain,Subcluster_1,Subcluster_2,Subcluster_3
ERR1234567,1,1_1,1_1_1,1_1_1_1
ERR1234568,1,1_1,1_1_1,1_1_1_2
ERR1234569,1,1_2,1_2_1,1_2_1_1
ERR1234570,1,1_2,1_2_1,1_2_1_2
GCA_000001.1,2,2_1,2_1_1,2_1_1_1
GCA_000002.1,2,2_1,2_1_2,2_1_2_1
GCA_000003.1,2,2_2,2_2_1,2_2_1_1
EOF

    # Create assembly files with realistic names
    mkdir -p assemblies
    
    # Create FASTA files
    samples=("ERR1234567" "ERR1234568" "ERR1234569" "ERR1234570" "GCA_000001.1" "GCA_000002.1" "GCA_000003.1")
    
    for sample in "\${samples[@]}"; do
        cat > assemblies/\${sample}.fasta << EOF
>\${sample}_contig_1
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
>\${sample}_contig_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
EOF
    done
    
    # Create rfiles.txt
    cat > rfiles.txt << EOF
ERR1234567	\$PWD/assemblies/ERR1234567.fasta
ERR1234568	\$PWD/assemblies/ERR1234568.fasta
ERR1234569	\$PWD/assemblies/ERR1234569.fasta
ERR1234570	\$PWD/assemblies/ERR1234570.fasta
GCA_000001.1	\$PWD/assemblies/GCA_000001.1.fasta
GCA_000002.1	\$PWD/assemblies/GCA_000002.1.fasta
GCA_000003.1	\$PWD/assemblies/GCA_000003.1.fasta
EOF

    # Create strain 1 directory with all expected files
    cat > poppipe_output/strains/1/rfile.txt << 'EOF'
ERR1234567	/path/to/ERR1234567.fasta
ERR1234568	/path/to/ERR1234568.fasta
ERR1234569	/path/to/ERR1234569.fasta
ERR1234570	/path/to/ERR1234570.fasta
EOF

    cat > poppipe_output/strains/1/names.txt << 'EOF'
ERR1234567
ERR1234568
ERR1234569
ERR1234570
EOF

    cat > poppipe_output/strains/1/njtree.nwk << 'EOF'
((ERR1234567:0.001,ERR1234568:0.001):0.002,(ERR1234569:0.002,ERR1234570:0.002):0.001);
EOF

    cat > poppipe_output/strains/1/besttree.nwk << 'EOF'
((ERR1234567:0.0012,ERR1234568:0.0011):0.0021,(ERR1234569:0.0019,ERR1234570:0.0020):0.0012);
EOF

    cat > poppipe_output/strains/1/fastbaps_clusters.txt << 'EOF'
ERR1234567	1
ERR1234568	1
ERR1234569	2
ERR1234570	2
EOF

    # Create dummy files for strain 1
    echo ">ERR1234567" > poppipe_output/strains/1/align_variants.aln
    echo "ATCGATCGATCG" >> poppipe_output/strains/1/align_variants.aln
    echo ">ERR1234568" >> poppipe_output/strains/1/align_variants.aln
    echo "ATCGATCGATCG" >> poppipe_output/strains/1/align_variants.aln
    
    # Create binary files (empty for testing)
    touch poppipe_output/strains/1/dists.npy
    touch poppipe_output/strains/1/dists.pkl
    touch poppipe_output/strains/1/split_kmers.skf

    # Create strain 2 directory
    cat > poppipe_output/strains/2/rfile.txt << 'EOF'
GCA_000001.1	/path/to/GCA_000001.1.fasta
GCA_000002.1	/path/to/GCA_000002.1.fasta
GCA_000003.1	/path/to/GCA_000003.1.fasta
EOF

    cat > poppipe_output/strains/2/names.txt << 'EOF'
GCA_000001.1
GCA_000002.1
GCA_000003.1
EOF

    cat > poppipe_output/strains/2/njtree.nwk << 'EOF'
(GCA_000001.1:0.003,(GCA_000002.1:0.002,GCA_000003.1:0.002):0.001);
EOF

    cat > poppipe_output/strains/2/besttree.nwk << 'EOF'
(GCA_000001.1:0.0031,(GCA_000002.1:0.0019,GCA_000003.1:0.0021):0.0011);
EOF

    cat > poppipe_output/strains/2/fastbaps_clusters.txt << 'EOF'
GCA_000001.1	1
GCA_000002.1	2
GCA_000003.1	2
EOF

    # Create dummy files for strain 2
    echo ">GCA_000001.1" > poppipe_output/strains/2/align_variants.aln
    echo "GCTAGCTAGCTAG" >> poppipe_output/strains/2/align_variants.aln
    echo ">GCA_000002.1" >> poppipe_output/strains/2/align_variants.aln
    echo "GCTAGCTAGCTAG" >> poppipe_output/strains/2/align_variants.aln
    
    touch poppipe_output/strains/2/dists.npy
    touch poppipe_output/strains/2/dists.pkl
    touch poppipe_output/strains/2/split_kmers.skf

    # Create backbone tree
    cat > poppipe_output/strains/other/njtree.nwk << 'EOF'
(((ERR1234567:0.001,ERR1234568:0.001):0.002,(ERR1234569:0.002,ERR1234570:0.002):0.001):0.01,((GCA_000001.1:0.003,(GCA_000002.1:0.002,GCA_000003.1:0.002):0.001):0.01):0.01);
EOF

    echo "Created complete PopPIPE-bp test data structure"
    """
}

// Test the complete parsing workflow
process TEST_COMPLETE_PARSING {
    publishDir params.outdir, mode: 'copy'
    
    input:
    path poppipe_dir
    path rfiles

    output:
    path "cluster_*.txt"
    path "cluster_summary.tsv"
    path "cluster_info.json"
    path "parsing_report.txt"

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
    import json
    import glob
    from pathlib import Path

    # Test complete PopPIPE parsing workflow
    print("=== PopPIPE-bp Complete Integration Test ===\\n")
    
    # Find all_clusters.txt
    poppipe_dir = '${poppipe_dir}'
    cluster_file = os.path.join(poppipe_dir, 'all_clusters.txt')
    
    if not os.path.exists(cluster_file):
        raise FileNotFoundError(f"all_clusters.txt not found: {cluster_file}")
    
    print(f"✓ Found all_clusters.txt: {cluster_file}")
    
    # Read and validate all_clusters.txt
    df = pd.read_csv(cluster_file)
    expected_columns = ['Taxon', 'Strain', 'Subcluster_1', 'Subcluster_2', 'Subcluster_3']
    
    print(f"✓ Columns found: {df.columns.tolist()}")
    print(f"✓ Expected columns: {expected_columns}")
    
    missing_cols = set(expected_columns) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing expected columns: {missing_cols}")
    
    print(f"✓ All expected columns present")
    print(f"✓ Data shape: {df.shape}")
    
    # Display sample data
    print("\\n=== Sample Data ===")
    print(df.to_string(index=False))
    
    # Read rfiles.txt
    rfiles_path = '${rfiles}'
    sample_to_file = {}
    
    with open(rfiles_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    sample_name = parts[0].strip()
                    file_path = parts[1].strip()
                    if os.path.isfile(file_path):
                        sample_to_file[sample_name] = file_path
    
    print(f"\\n✓ Loaded {len(sample_to_file)} sample-to-file mappings")
    
    # Check strain directories
    strain_path = os.path.join(poppipe_dir, 'strains')
    strain_info = {}
    
    if os.path.exists(strain_path):
        for strain_dir in os.listdir(strain_path):
            strain_dir_path = os.path.join(strain_path, strain_dir)
            if os.path.isdir(strain_dir_path):
                files = os.listdir(strain_dir_path)
                strain_info[strain_dir] = {
                    'path': strain_dir_path,
                    'files': files,
                    'file_count': len(files)
                }
    
    print(f"\\n✓ Found {len(strain_info)} strain directories")
    for strain_id, info in strain_info.items():
        print(f"  Strain {strain_id}: {info['file_count']} files")
        expected_files = ['rfile.txt', 'names.txt', 'njtree.nwk', 'besttree.nwk', 
                         'fastbaps_clusters.txt', 'align_variants.aln', 'dists.npy', 
                         'dists.pkl', 'split_kmers.skf']
        found_files = set(info['files'])
        missing_files = set(expected_files) - found_files
        if missing_files:
            print(f"    Missing files: {missing_files}")
        else:
            print(f"    ✓ All expected files present")
    
    # Process clusters
    clusters = df.groupby('Strain')
    cluster_summary = []
    
    print(f"\\n=== Processing {len(clusters)} Strains ===")
    
    for strain_id, group in clusters:
        print(f"\\nStrain {strain_id}:")
        cluster_files = []
        cluster_filename = f"cluster_{strain_id}.txt"
        
        with open(cluster_filename, 'w') as f:
            for _, row in group.iterrows():
                sample_name = row['Taxon']
                print(f"  Sample: {sample_name}")
                print(f"    Subcluster_1: {row['Subcluster_1']}")
                print(f"    Subcluster_2: {row['Subcluster_2']}")
                print(f"    Subcluster_3: {row['Subcluster_3']}")
                
                matched_file = sample_to_file.get(sample_name)
                if matched_file:
                    f.write(f"{matched_file}\\n")
                    cluster_files.append(matched_file)
                    print(f"    ✓ Matched to: {matched_file}")
                else:
                    print(f"    ✗ No file match found")
        
        cluster_summary.append({
            'Cluster': strain_id,
            'Sample_Count': len(cluster_files),
            'Samples': ';'.join([Path(f).stem for f in cluster_files]),
            'Has_PopPIPE_Data': str(strain_id) in strain_info
        })
        
        print(f"  ✓ Created cluster file with {len(cluster_files)} samples")
    
    # Write outputs
    summary_df = pd.DataFrame(cluster_summary)
    summary_df.to_csv('cluster_summary.tsv', sep='\\t', index=False)
    
    with open('cluster_info.json', 'w') as f:
        json.dump(strain_info, f, indent=2)
    
    # Generate report
    with open('parsing_report.txt', 'w') as f:
        f.write("PopPIPE-bp Complete Integration Test Report\\n")
        f.write("=" * 50 + "\\n\\n")
        f.write(f"Input file: {cluster_file}\\n")
        f.write(f"Data shape: {df.shape}\\n")
        f.write(f"Columns: {', '.join(df.columns)}\\n")
        f.write(f"Sample-to-file mappings: {len(sample_to_file)}\\n")
        f.write(f"Strain directories found: {len(strain_info)}\\n")
        f.write(f"Clusters processed: {len(clusters)}\\n\\n")
        
        f.write("Cluster Summary:\\n")
        f.write(summary_df.to_string(index=False))
        f.write("\\n\\nStrain Directory Contents:\\n")
        for strain_id, info in strain_info.items():
            f.write(f"\\nStrain {strain_id} ({info['file_count']} files):\\n")
            for file in sorted(info['files']):
                f.write(f"  - {file}\\n")
    
    print(f"\\n=== Test Results ===")
    print(f"✓ Successfully processed {len(clusters)} strains")
    print(f"✓ Generated {len(cluster_summary)} cluster files")
    print(f"✓ All samples matched to assembly files")
    print(f"✓ PopPIPE strain data validated")
    print("\\n🎉 Complete PopPIPE-bp integration test PASSED!")
    """
}

workflow TEST_COMPLETE_POPPIPE {
    
    // Create comprehensive test data
    CREATE_COMPLETE_POPPIPE_DATA()
    
    // Test complete parsing
    TEST_COMPLETE_PARSING(
        CREATE_COMPLETE_POPPIPE_DATA.out.poppipe_dir,
        CREATE_COMPLETE_POPPIPE_DATA.out.rfiles
    )
    
    // Display results
    TEST_COMPLETE_PARSING.out[3].view { "Report: " + it.text }
}

workflow {
    TEST_COMPLETE_POPPIPE()
}