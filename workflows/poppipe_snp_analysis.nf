/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    POPPIPE SNP ANALYSIS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Integration workflow for PopPIPE-bp outputs with SNP analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSNPS.initialise(params, log)

// Check input path parameters
def checkPathParamList = [ params.poppipe_output, params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.poppipe_output) { ch_poppipe_output = file(params.poppipe_output) } else { exit 1, 'PopPIPE output directory not specified!' }
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input assembly directory not specified!' }
if (params.ref) { ch_ref_input = file(params.ref) } else { ch_ref_input = [] }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Local modules
//
include { PARSE_POPPIPE_OUTPUT                             } from "../modules/local/parse_poppipe_output/main"
include { TREE_GRAFT                                       } from "../modules/local/tree_graft/main"
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

workflow POPPIPE_SNP_ANALYSIS {

    // SETUP: Define empty channels to concatenate certain outputs
    ch_versions             = Channel.empty()
    ch_qc_filecheck         = Channel.empty()
    ch_output_summary_files = Channel.empty()

    /*
    ================================================================================
                            Preprocess input data
    ================================================================================
    */

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
                        Parse PopPIPE-bp outputs
    ================================================================================
    */

    log.info "Parsing PopPIPE-bp output directory: ${params.poppipe_output}"
    log.info "Input source: ${params.input}"
    
    // Parse PopPIPE output and create cluster files
    // The input can be either a TSV file (rfiles.txt) or a directory
    PARSE_POPPIPE_OUTPUT (
        ch_poppipe_output,
        ch_input,
        Channel.value(params.input)  // Pass original input path for relative path resolution
    )
    ch_versions = ch_versions.mix(PARSE_POPPIPE_OUTPUT.out.versions)

    /*
    ================================================================================
                        Run SNP analysis on each cluster
    ================================================================================
    */

    // Create channels for each cluster
    PARSE_POPPIPE_OUTPUT.out.cluster_files
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

    /*
    ================================================================================
                        Tree grafting and integration
    ================================================================================
    */

    // Collect all phylogenetic trees from clusters
    ch_cluster_trees = CLUSTER_SNP_ANALYSIS.out.phylogenetic_tree
        .map { meta, tree -> tree }
        .collect()

    // Perform tree grafting using PopPIPE algorithm
    TREE_GRAFT (
        ch_cluster_trees,
        ch_poppipe_output,
        PARSE_POPPIPE_OUTPUT.out.cluster_info
    )
    ch_versions = ch_versions.mix(TREE_GRAFT.out.versions)

    /*
    ================================================================================
                        Collect outputs and summaries
    ================================================================================
    */

    // Collect outputs from all clusters
    ch_output_summary_files = ch_output_summary_files.mix(
        CLUSTER_SNP_ANALYSIS.out.distance_pairs.map { meta, file -> file },
        CLUSTER_SNP_ANALYSIS.out.distance_matrix.map { meta, file -> file },
        CLUSTER_SNP_ANALYSIS.out.masked_distance_matrix.map { meta, file -> file },
        PARSE_POPPIPE_OUTPUT.out.summary,
        TREE_GRAFT.out.tree,
        TREE_GRAFT.out.log
    )

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
    grafted_tree = TREE_GRAFT.out.tree
    cluster_summaries = PARSE_POPPIPE_OUTPUT.out.summary
    cluster_trees = CLUSTER_SNP_ANALYSIS.out.phylogenetic_tree
    cluster_alignments = CLUSTER_SNP_ANALYSIS.out.core_alignment
    distance_matrices = CLUSTER_SNP_ANALYSIS.out.distance_matrix
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