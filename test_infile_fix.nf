#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INFILE_HANDLING_UNIX } from './modules/local/infile_handling_unix/main'

workflow {
    // Create a test FASTA file
    Channel.of("test_sample")
        .map { sample_id ->
            def meta = [:]
            meta.id = sample_id
            
            // Create a simple test FASTA content
            def fasta_content = ">test_sequence\nATCGATCGATCGATCG\n"
            def test_file = file("${sample_id}.fasta")
            test_file.text = fasta_content
            
            return [meta, test_file]
        }
        .set { ch_test_input }
    
    INFILE_HANDLING_UNIX(ch_test_input)
    
    INFILE_HANDLING_UNIX.out.qc_filecheck.view { "QC Check: $it" }
    INFILE_HANDLING_UNIX.out.input_files.view { "Input Files: $it" }
}