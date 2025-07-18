/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSNPS.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet or directory not specified!' }
if (params.ref) { ch_ref_input = file(params.ref) } else { ch_ref_input = [] }

// Check for PopPUNK cluster file
if (params.poppunk_clusters) { 
    ch_poppunk_clusters = file(params.poppunk_clusters) 
    if (!ch_poppunk_clusters.exists()) {
        exit 1, "PopPUNK cluster file not found: ${params.poppunk_clusters}"
    }
} else { 
    ch_poppunk_clusters = [] 
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { PARSE_POPPUNK_CLUSTERS                           } from "../modules/local/parse_poppunk_clusters/main"
include { INFILE_HANDLING_UNIX                             } from "../modules/local/infile_handling_unix/main"
include { INFILE_HANDLING_UNIX as REF_INFILE_HANDLING_UNIX } from "../modules/local/infile_handling_unix/main"
include { CONVERT_TSV_TO_EXCEL_PYTHON                      } from "../modules/local/convert_tsv_to_excel_python/main"
include { CREATE_EXCEL_RUN_SUMMARY_PYTHON                  } from "../modules/local/create_excel_run_summary_python/main"

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                                      } from "../subworkflows/local/input_check"
include { INPUT_CHECK as REF_INPUT_CHECK                   } from "../subworkflows/local/input_check"
include { CLUSTER_SNP_ANALYSIS                             } from "../subworkflows/local/cluster_snp_analysis"
include { ASSEMBLY_SNPS as REGULAR_ASSEMBLY_SNPS           } from "./assembly_snps"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE CHANNELS FOR INPUT PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if ( toLower(params.snp_package) == "parsnp" ) {
    ch_snp_package = "Parsnp"
} else {
    ch_snp_package = "Parsnp"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Convert input to lowercase
def toLower(it) {
    it.toString().toLowerCase()
}

// Check QC filechecks for a failure
def qcfilecheck(process, qcfile, inputfile) {
    qcfile.map{ meta, file -> [ meta, [file] ] }
            .join(inputfile)
            .map{ meta, qc, input ->
                data = []
                qc.flatten().each{ data += it.readLines() }

                if ( data.any{ it.contains('FAIL') } ) {
                    line = data.last().split('\t')
                    if (line.first() != "NaN") {
                        log.warn("${line[1]} QC check failed during process ${process} for sample ${line.first()}")
                    } else {
                        log.warn("${line[1]} QC check failed during process ${process}")
                    }
                } else {
                    [ meta, input ]
                }
            }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLY_SNPS_CLUSTERED {

    // SETUP: Define empty channels to concatenate certain outputs
    ch_versions             = Channel.empty()
    ch_qc_filecheck         = Channel.empty()
    ch_output_summary_files = Channel.empty()

    /*
    ================================================================================
                            Preprocess input data
    ================================================================================
    */

    // SUBWORKFLOW: Check input for samplesheet or pull inputs from directory
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Handle reference file if provided
    if (params.ref) {
        REF_INPUT_CHECK (
            ch_ref_input
        )
        ch_versions = ch_versions.mix(REF_INPUT_CHECK.out.versions)

        REF_INFILE_HANDLING_UNIX (
            REF_INPUT_CHECK.out.input_files
        )
        ch_versions = ch_versions.mix(REF_INFILE_HANDLING_UNIX.out.versions)
        ch_qc_filecheck = ch_qc_filecheck.concat(REF_INFILE_HANDLING_UNIX.out.qc_filecheck)

        ch_reference_file = REF_INFILE_HANDLING_UNIX.out.input_files
                                .map { meta, file -> file }
                                .first()
    } else {
        ch_reference_file = Channel.value(null)
    }

    /*
    ================================================================================
                        Choose processing mode: Clustered vs Regular
    ================================================================================
    */

    if (params.poppunk_clusters) {
        log.info "Running in PopPUNK cluster mode"
        
        // Parse PopPUNK clusters and create file lists for each cluster
        // Create a channel with file paths instead of staging files
        ch_input_file_paths = INPUT_CHECK.out.input_files
            .map { meta, file -> file.toString() }
            .collect()
            .map { paths -> paths.join('\n') }
        
        PARSE_POPPUNK_CLUSTERS (
            ch_poppunk_clusters,
            ch_input_file_paths
        )
        ch_versions = ch_versions.mix(PARSE_POPPUNK_CLUSTERS.out.versions)

        // Read cluster files and create channels for each cluster
        PARSE_POPPUNK_CLUSTERS.out.cluster_files
            .flatten()
            .map { cluster_file ->
                def cluster_id = cluster_file.baseName.replaceAll('cluster_', '')
                def file_lines = cluster_file.readLines()
                def files = file_lines.findAll { it.trim() != "" }.collect { file(it.trim()) }
                [ cluster_id, files ]
            }
            .filter { cluster_id, files -> files.size() >= params.min_cluster_size }
            .set { ch_clusters }

        // Run SNP analysis on each cluster
        CLUSTER_SNP_ANALYSIS (
            ch_clusters,
            ch_reference_file,
            Channel.value(ch_snp_package)
        )
        ch_versions = ch_versions.mix(CLUSTER_SNP_ANALYSIS.out.versions)
        ch_qc_filecheck = ch_qc_filecheck.concat(CLUSTER_SNP_ANALYSIS.out.qc_filecheck)
        
        // Collect outputs from all clusters
        ch_output_summary_files = ch_output_summary_files.mix(
            CLUSTER_SNP_ANALYSIS.out.distance_pairs.map { meta, file -> file },
            CLUSTER_SNP_ANALYSIS.out.distance_matrix.map { meta, file -> file },
            CLUSTER_SNP_ANALYSIS.out.masked_distance_matrix.map { meta, file -> file }
        )

    } else {
        log.info "Running in regular mode (all samples together)"
        
        // Run the original workflow on all samples
        REGULAR_ASSEMBLY_SNPS()
        ch_versions = ch_versions.mix(REGULAR_ASSEMBLY_SNPS.out.versions)
    }

    /*
    ================================================================================
                        Collect QC information
    ================================================================================
    */

    // Collect QC file check information
    ch_qc_filecheck = ch_qc_filecheck
                        .map{ meta, file -> file }
                        .collectFile(
                            name:       "Summary.QC_File_Checks.tsv",
                            keepHeader: true,
                            storeDir:   "${params.outdir}/Summaries",
                            sort:       'index'
                        )

    ch_output_summary_files = ch_output_summary_files.mix(ch_qc_filecheck.collect())

    /*
    ================================================================================
                        Convert TSV outputs to Excel XLSX
    ================================================================================
    */

    if (params.create_excel_outputs) {
        CREATE_EXCEL_RUN_SUMMARY_PYTHON (
            ch_output_summary_files.collect()
        )
        ch_versions = ch_versions.mix(CREATE_EXCEL_RUN_SUMMARY_PYTHON.out.versions)

        CONVERT_TSV_TO_EXCEL_PYTHON (
            CREATE_EXCEL_RUN_SUMMARY_PYTHON.out.summary
        )
        ch_versions = ch_versions.mix(CONVERT_TSV_TO_EXCEL_PYTHON.out.versions)
    }

    /*
    ================================================================================
                        Collect version information
    ================================================================================
    */

    // Collect version information
    ch_versions
        .unique()
        .collectFile(
            name:     "software_versions.yml",
            storeDir: params.tracedir
        )

    emit:
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/