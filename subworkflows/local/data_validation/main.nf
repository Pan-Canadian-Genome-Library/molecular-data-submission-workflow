// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { VALIDATION_FILEINTEGRITY } from '../../../modules/local/validation/fileintegrity/main'
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

    VALIDATION_FILEINTEGRITY ( ch_integrity_input )
    ch_versions = ch_versions.mix(VALIDATION_FILEINTEGRITY.out.versions.first())

    // Reconstruct validated_payload_files channel with final meta.status
    // Get final results from file integrity validation with updated meta.status
    ch_validated_payload_files = VALIDATION_FILEINTEGRITY.out.ch_payload_files
        .join(VALIDATION_FILEINTEGRITY.out.status, by: 0)
        .map { meta, payload, files, status_file ->
            // Update meta with final validation status
            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()
            updated_meta.status = status_value
            [updated_meta, payload, files]
        }

    // Collect all status files from all validation steps
    ch_all_status = Channel.empty()
        .mix(VALIDATION_METADATA.out.status)
        .mix(VALIDATION_CROSSCHECK.out.status)
        .mix(VALIDATION_FILEINTEGRITY.out.status)

    emit:
    // Output channels as specified
    validated_payload_files = ch_validated_payload_files  // channel: [ val(meta), payload(json), [files] ]
    status                  = ch_all_status               // channel: [ val(meta), *_status.yml ] - mixed from all validation steps

    versions = ch_versions                                // channel: [ versions.yml ]
}
