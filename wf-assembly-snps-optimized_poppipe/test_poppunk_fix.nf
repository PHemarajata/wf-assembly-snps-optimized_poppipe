#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Test script to validate the PARSE_POPPUNK_CLUSTERS fix
include { INPUT_CHECK } from "./subworkflows/local/input_check"
include { PARSE_POPPUNK_CLUSTERS } from "./modules/local/parse_poppunk_clusters/main"

workflow {
    // Create test data - simulate the file collision scenario
    // This would normally come from your actual input
    
    // For testing, let's create a simple channel structure
    ch_test_files = Channel.of(
        [['id': 'sample1'], file('/path/to/GCA_963562795.1.fasta')],
        [['id': 'sample2'], file('/path/to/different/GCA_963562795.1.fasta')],
        [['id': 'sample3'], file('/path/to/IE-0014_S8_L001-SPAdes.fasta')]
    )
    
    // Create test PopPUNK assignments file
    ch_poppunk = Channel.of(file('test_poppunk_clusters.csv'))
    
    // Test the file path collection
    ch_input_file_paths = ch_test_files
        .map { meta, file -> file.toString() }
        .collect()
        .map { paths -> paths.join('\n') }
        .view { "File paths: ${it}" }
    
    // This would be the actual call (commented out for testing)
    // PARSE_POPPUNK_CLUSTERS (
    //     ch_poppunk,
    //     ch_input_file_paths
    // )
}