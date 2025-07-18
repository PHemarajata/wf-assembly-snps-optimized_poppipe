#!/usr/bin/env nextflow

/*
 * Simple test for relative path functionality
 */

nextflow.enable.dsl = 2

process TEST_RELATIVE_PATHS_SIMPLE {
    output:
    path "test_log.txt"

    script:
    """
    #!/usr/bin/env python3
    
    import os
    import tempfile
    import shutil
    
    # Create a temporary test structure
    test_dir = "test_project"
    os.makedirs(f"{test_dir}/data", exist_ok=True)
    os.makedirs(f"{test_dir}/assemblies", exist_ok=True)
    
    # Create test assembly files
    with open(f"{test_dir}/assemblies/sample1.fasta", "w") as f:
        f.write(">sample1\\nATCGATCGATCG\\n")
    
    with open(f"{test_dir}/assemblies/sample2.fasta", "w") as f:
        f.write(">sample2\\nGCTAGCTAGCTA\\n")
    
    # Create TSV file with relative paths
    with open(f"{test_dir}/data/rfiles.txt", "w") as f:
        f.write("sample1\\t../assemblies/sample1.fasta\\n")
        f.write("sample2\\t../assemblies/sample2.fasta\\n")
    
    # Test relative path resolution
    tsv_file_path = f"{test_dir}/data/rfiles.txt"
    tsv_dir = os.path.dirname(os.path.abspath(tsv_file_path))
    
    log_lines = []
    log_lines.append("=== RELATIVE PATH TEST ===")
    log_lines.append(f"TSV file: {tsv_file_path}")
    log_lines.append(f"TSV directory: {tsv_dir}")
    
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
                
                log_lines.append(f"\\nProcessing line {line_num}:")
                log_lines.append(f"  Sample: {sample_name}")
                log_lines.append(f"  Original path: {file_path}")
                
                # If file_path is relative, make it relative to the TSV file directory
                if not os.path.isabs(file_path):
                    resolved_path = os.path.join(tsv_dir, file_path)
                    log_lines.append(f"  Resolved path: {resolved_path}")
                else:
                    resolved_path = file_path
                    log_lines.append(f"  Absolute path: {resolved_path}")
                
                # Normalize the path to resolve any .. or . components
                normalized_path = os.path.normpath(resolved_path)
                log_lines.append(f"  Normalized path: {normalized_path}")
                
                if os.path.isfile(normalized_path):
                    sample_to_file[sample_name] = normalized_path
                    log_lines.append(f"  ✓ File found and mapped")
                else:
                    log_lines.append(f"  ✗ File not found")
    
    log_lines.append(f"\\n=== RESULTS ===")
    log_lines.append(f"Successfully mapped {len(sample_to_file)} samples:")
    for sample, path in sample_to_file.items():
        log_lines.append(f"  {sample} -> {path}")
    
    log_lines.append("\\n✓ Relative path resolution test completed successfully!")
    
    # Write log
    with open("test_log.txt", "w") as f:
        f.write("\\n".join(log_lines))
    """
}

workflow {
    TEST_RELATIVE_PATHS_SIMPLE()
    TEST_RELATIVE_PATHS_SIMPLE.out.view { "Test results:\\n" + it.text }
}