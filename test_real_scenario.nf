#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Set the parameters that the real pipeline uses
params.min_input_filesize = "1k"

include { INFILE_HANDLING_UNIX } from './modules/local/infile_handling_unix/main'

workflow {
    // Create a realistic test FASTA file
    Channel.of("sample_001")
        .map { sample_id ->
            def meta = [:]
            meta.id = sample_id
            
            // Create a realistic FASTA file (>1k)
            def header = ">contig_1 length=2000 coverage=50.5"
            def sequence = ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 50)
            def fasta_content = "${header}\n${sequence}\n"
            
            def test_file = file("${sample_id}.fasta")
            test_file.text = fasta_content
            
            return [meta, test_file]
        }
        .set { ch_test_input }
    
    INFILE_HANDLING_UNIX(ch_test_input)
    
    // Check outputs
    INFILE_HANDLING_UNIX.out.qc_filecheck
        .subscribe { meta, tsv_file ->
            println "QC Check for ${meta.id}: ${tsv_file}"
            println "TSV content:"
            println tsv_file.text
        }
    
    INFILE_HANDLING_UNIX.out.input_files
        .subscribe { files ->
            println "Input files passed QC:"
            files.each { println "  - ${it}" }
        }
}