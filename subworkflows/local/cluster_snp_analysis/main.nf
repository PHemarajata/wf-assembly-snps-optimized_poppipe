/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLUSTER SNP ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INFILE_HANDLING_UNIX                             } from "../../../modules/local/infile_handling_unix/main"
include { CORE_GENOME_ALIGNMENT_PARSNP                     } from "../../../modules/local/core_genome_alignment_parsnp/main"
include { CONVERT_GINGR_TO_FASTA_HARVESTTOOLS              } from "../../../modules/local/convert_gingr_to_fasta_harvesttools/main"
include { CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS           } from "../../../modules/local/calculate_pairwise_distances_snp_dists/main"
include { CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS             } from "../../../modules/local/create_snp_distance_matrix_snp_dists/main"
include { MASK_RECOMBINANT_POSITIONS_BIOPYTHON             } from "../../../modules/local/mask_recombinant_positions_biopython/main"
include { CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS      } from "../../../modules/local/create_masked_snp_distance_matrix_snp_dists/main"
include { BUILD_PHYLOGENETIC_TREE_PARSNP                   } from "../../../modules/local/build_phylogenetic_tree_parsnp/main"
include { RECOMBINATION                                    } from "../recombination"

workflow CLUSTER_SNP_ANALYSIS {
    take:
    cluster_files    // channel: [ cluster_id, [files] ]
    reference_file   // val: reference file (optional)
    snp_package      // val: snp package name

    main:
    ch_versions = Channel.empty()
    ch_qc_filecheck = Channel.empty()

    // Process each cluster
    cluster_files
        .combine(snp_package)
        .map { cluster_id, files, snp_pkg ->
            def meta = [:]
            meta['id'] = "cluster_${cluster_id}"
            meta['cluster_id'] = cluster_id
            meta['snp_package'] = snp_pkg
            [ meta, files ]
        }
        .set { ch_cluster_input }

    // Check input files meet size criteria for each cluster
    INFILE_HANDLING_UNIX (
        ch_cluster_input.flatMap { meta, files ->
            files.collect { file ->
                def file_meta = [:]
                file_meta.putAll(meta)
                file_meta['id'] = "${meta.cluster_id}_${file.baseName}"
                [ file_meta, file ]
            }
        }
    )
    ch_versions = ch_versions.mix(INFILE_HANDLING_UNIX.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(INFILE_HANDLING_UNIX.out.qc_filecheck)

    // Group files back by cluster
    INFILE_HANDLING_UNIX.out.input_files
        .map { meta, file ->
            def cluster_id = meta.cluster_id
            def snp_pkg = meta.snp_package
            [ cluster_id, file, snp_pkg ]
        }
        .groupTuple()
        .map { cluster_id, files, snp_pkgs ->
            def meta = [:]
            meta['id'] = "cluster_${cluster_id}"
            meta['cluster_id'] = cluster_id
            meta['snp_package'] = snp_pkgs[0] // All should be the same
            [ meta, files ]
        }
        .set { ch_cluster_files }

    // Handle reference file for each cluster
    ch_cluster_files
        .combine(reference_file)
        .map { meta, files, ref_file ->
            if (ref_file != null) {
                [ meta, files, ref_file ]
            } else {
                // Use largest file in each cluster as reference
                def sorted_files = files.sort { it.size() }.reverse()
                def ref_file_local = sorted_files[0]
                def remaining_files = sorted_files.drop(1)
                [ meta, remaining_files, ref_file_local ]
            }
        }
        .set { ch_cluster_with_ref }

    // Run ParSNP for each cluster
    CORE_GENOME_ALIGNMENT_PARSNP (
        ch_cluster_with_ref.map { meta, files, ref -> [ meta, files ] },
        ch_cluster_with_ref.map { meta, files, ref -> [ meta, ref ] }
    )
    ch_versions = ch_versions.mix(CORE_GENOME_ALIGNMENT_PARSNP.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(CORE_GENOME_ALIGNMENT_PARSNP.out.qc_filecheck)

    // Convert Gingr to FastA
    CONVERT_GINGR_TO_FASTA_HARVESTTOOLS (
        CORE_GENOME_ALIGNMENT_PARSNP.out.output
    )
    ch_versions = ch_versions.mix(CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.qc_filecheck)

    // Calculate distances
    CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS (
        CORE_GENOME_ALIGNMENT_PARSNP.out.output
    )
    ch_versions = ch_versions.mix(CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS.out.versions)

    CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS (
        CORE_GENOME_ALIGNMENT_PARSNP.out.output
    )
    ch_versions = ch_versions.mix(CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS.out.versions)

    // Recombination analysis
    RECOMBINATION (
        CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.core_alignment,
        CORE_GENOME_ALIGNMENT_PARSNP.out.output
    )
    ch_versions = ch_versions.mix(RECOMBINATION.out.versions)

    // Mask recombinant positions
    MASK_RECOMBINANT_POSITIONS_BIOPYTHON (
        RECOMBINATION.out.recombinants,
        CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.core_alignment.collect()
    )
    ch_versions = ch_versions.mix(MASK_RECOMBINANT_POSITIONS_BIOPYTHON.out.versions)

    // Create masked distance matrix
    CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS (
        MASK_RECOMBINANT_POSITIONS_BIOPYTHON.out.masked_alignment
    )
    ch_versions = ch_versions.mix(CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS.out.versions)

    // Build phylogenetic tree
    BUILD_PHYLOGENETIC_TREE_PARSNP (
        MASK_RECOMBINANT_POSITIONS_BIOPYTHON.out.masked_alignment
    )
    ch_versions = ch_versions.mix(BUILD_PHYLOGENETIC_TREE_PARSNP.out.versions)
    ch_qc_filecheck = ch_qc_filecheck.concat(BUILD_PHYLOGENETIC_TREE_PARSNP.out.qc_filecheck)

    emit:
    versions = ch_versions
    qc_filecheck = ch_qc_filecheck
    distance_pairs = CALCULATE_PAIRWISE_DISTANCES_SNP_DISTS.out.snp_distances
    distance_matrix = CREATE_SNP_DISTANCE_MATRIX_SNP_DISTS.out.distance_matrix
    masked_distance_matrix = CREATE_MASKED_SNP_DISTANCE_MATRIX_SNP_DISTS.out.distance_matrix
    phylogenetic_tree = BUILD_PHYLOGENETIC_TREE_PARSNP.out.tree
    core_alignment = CONVERT_GINGR_TO_FASTA_HARVESTTOOLS.out.core_alignment
    masked_alignment = MASK_RECOMBINANT_POSITIONS_BIOPYTHON.out.masked_alignment
}