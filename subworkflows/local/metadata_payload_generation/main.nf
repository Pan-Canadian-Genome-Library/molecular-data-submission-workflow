// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { PAYLOAD_GENERATE   } from '../../../modules/local/payload/generate/main'
include { PAYLOAD_VALIDATE   } from '../../../modules/local/payload/validate/main'

workflow METADATA_PAYLOAD_GENERATION {

    take:
    molecular_metadata_files // channel: [val(meta.id=submitter_analysis_id, meta.type=analysisType, meta.study=studyId), path(file_meta.tsv), path(analysis_meta.tsv), path(workflow_meta.tsv), [data_files]] per analysis

    main:

    ch_versions = Channel.empty()

    // Step 1: Generate payloads (always runs)
    PAYLOAD_GENERATE ( molecular_metadata_files )
    ch_versions = ch_versions.mix(PAYLOAD_GENERATE.out.versions.first())

    // Step 2: Validate payloads - update meta.status based on generation result
    ch_validate_input = PAYLOAD_GENERATE.out.payload_files
        .join(PAYLOAD_GENERATE.out.status, by: 0)
        .map { meta, payload_file, data_files, status_file ->
            // Read status from YAML file and update meta
            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()
            updated_meta.status = status_value
            [updated_meta, payload_file, data_files]
        }

    PAYLOAD_VALIDATE ( ch_validate_input )
    ch_versions = ch_versions.mix(PAYLOAD_VALIDATE.out.versions.first())

    // Get final results from validation with updated meta.status
    validated_payload_files = PAYLOAD_VALIDATE.out.validated_payload_files
        .join(PAYLOAD_VALIDATE.out.status, by: 0)
        .map { meta, payload_file, data_files, status_file ->
            // Read status from YAML file and update meta
            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()
            // Only successful if validation passed AND previous status was pass
            updated_meta.status = (status_value == 'pass' && meta.status == 'pass') ? 'pass' : 'failed'
            [updated_meta, payload_file, data_files]
        }

    // Simple filtering based on meta.status
    successful_analyses = validated_payload_files
        .filter { meta, _payload_file, _data_files -> 
            meta.status == 'pass' 
        }
    
    failed_analyses = validated_payload_files
        .filter { meta, _payload_file, _data_files ->
            meta.status == 'failed'
        }

    // Collect all status files from payload generation pipeline
    ch_all_status = Channel.empty()
    ch_all_status = ch_all_status.mix(PAYLOAD_GENERATE.out.status)
    ch_all_status = ch_all_status.mix(PAYLOAD_VALIDATE.out.status)

    emit:
    // Main outputs
    all_analyses     = validated_payload_files         // channel: [ val(meta_with_status), payload.json, [data_files] ]
    successful       = successful_analyses             // channel: [ val(meta), payload.json, [data_files] ] where meta.status == 'pass'
    failed           = failed_analyses                 // channel: [ val(meta), payload.json, [data_files] ] where meta.status == 'failed'
    
    // Status and versions for full workflow aggregation
    status           = ch_all_status                   // channel: [ val(meta), path(*_status.yml) ] - individual process status files
    versions         = ch_versions                     // channel: [ versions.yml ]
}
