// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { FILE_INTEGRITY         } from '../file_integrity/main'
include { VALIDATION_METADATA      } from '../../../modules/local/validation/metadata/main'
include { VALIDATION_CROSSCHECK    } from '../../../modules/local/validation/crosscheck/main'

workflow DATA_VALIDATION {

    take:
    // Input channels as specified
    ch_payload_files_biospecimen       // channel: [ val(meta), payload(json), [files], specimen(path), sample(path), experiment(path), read_group(path)] - from molecular_payload_generation

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow
    // Validation execution order: 1. metadata → 2. crosscheck → 3. fileintegrity

    // Step 1: Validate metadata consistency with analysis payload
    // Pass the two input channels directly to VALIDATION_METADATA
    VALIDATION_METADATA ( ch_payload_files_biospecimen )
    ch_versions = ch_versions.mix(VALIDATION_METADATA.out.versions.first())

    // Step 2: MD5 checksum cross-check validation
    // Update meta.status based on metadata validation result before passing to crosscheck
    ch_crosscheck_input = VALIDATION_METADATA.out.ch_payload_files
        .join(VALIDATION_METADATA.out.status, by: 0)
        .map { meta, payload, payload_files, status_file ->
            // Read status from YAML file and update meta
            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()
            updated_meta.status = status_value
            [updated_meta, payload, payload_files]
        }
    
    VALIDATION_CROSSCHECK ( ch_crosscheck_input )
    ch_versions = ch_versions.mix(VALIDATION_CROSSCHECK.out.versions.first())

    // Step 3: Validate file integrity (BAM, CRAM, VCF, FASTQ format checks)
    // Update meta.status based on crosscheck validation result before passing to file integrity
    ch_integrity_input = VALIDATION_CROSSCHECK.out.ch_payload_files
        .join(VALIDATION_CROSSCHECK.out.status, by: 0)
        .map { meta, payload, files, status_file ->
            // Read status from YAML file and update meta
            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()
            updated_meta.status = status_value
            // Keep the same structure as crosscheck output for fileintegrity
            [updated_meta, payload, files]
        }

    FILE_INTEGRITY ( ch_integrity_input )
    ch_versions = ch_versions.mix(FILE_INTEGRITY.out.versions.first())

    // Aggregate status from FILE_INTEGRITY validation which may contain multiple status files
    // Group status files by meta.id to aggregate all validation results per sample
    ch_file_integrity_aggregated_status = FILE_INTEGRITY.out.status
        .map { meta, status_file -> [meta.id, meta, status_file] }
        .groupTuple(by: 0)
        .map { _id, metas, status_files ->
            // Use first meta (should be identical for same sample)
            def meta = metas[0]
            
            // Aggregate status from all file validation results
            def overall_status = "pass"  // Default status
            
            // Check each status file for failures
            status_files.each { status_file ->
                if (status_file.exists()) {
                    def content = status_file.text.toLowerCase()
                    // Check for failure indicators in status file content
                    if (content.contains('status: "failed"')) {
                        overall_status = "failed"
                    }
                }
            }
            
            // Handle upstream failures - preserve existing failed status
            if (meta.status == "failed") {
                overall_status = "failed"
            }
            
            // Create updated meta with aggregated status
            def updated_meta = meta.clone()
            updated_meta.status = overall_status

            [updated_meta.id, updated_meta]
        }

    // Reconstruct validated_payload_files channel with final aggregated meta.status
    // Get final results from file integrity validation with updated meta.status
    ch_validated_payload_files = FILE_INTEGRITY.out.ch_payload_files
        .map { meta, payload, files -> [meta.id, meta, payload, files] }
        .join(ch_file_integrity_aggregated_status, by: 0)
        .map { _id, _original_meta, payload, files, updated_meta ->
            // Use the updated meta with aggregated status
            [updated_meta, payload, files]
        }

    // Collect all status files from all validation steps
    // Use original FILE_INTEGRITY status files since we only update meta.status internally
    ch_all_status = Channel.empty()
        .mix(VALIDATION_METADATA.out.status)
        .mix(VALIDATION_CROSSCHECK.out.status)
        .mix(FILE_INTEGRITY.out.status)

    // Debug output: View the final validated payload files channel
    ch_validated_payload_files.view { meta, payload, files ->
        "DATA_VALIDATION OUTPUT: meta.id=${meta.id}, meta.status=${meta.status}, payload=${payload.name}, files=[${files.collect{it.name}.join(', ')}] (${files.size()} total)"
    }

    emit:
    // Output channels as specified
    validated_payload_files = ch_validated_payload_files  // channel: [ val(meta), payload(json), [files] ]
    status                  = ch_all_status               // channel: [ val(meta), *_status.yml ] - mixed from all validation steps

    versions = ch_versions                                // channel: [ versions.yml ]
}
