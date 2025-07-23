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
    ch_status = Channel.empty()

    // Step 1: Generate payloads (always runs)
    PAYLOAD_GENERATE ( molecular_metadata_files )
    ch_versions = ch_versions.mix(PAYLOAD_GENERATE.out.versions.first())
    ch_status = ch_status.mix(PAYLOAD_GENERATE.out.status)

    // Update meta.status after generation using status.yml
    generation_with_status = PAYLOAD_GENERATE.out.payload_files
        .map { meta, payload_file, data_files -> [meta.id, meta, payload_file, data_files] }
        .join(
            PAYLOAD_GENERATE.out.status.map { meta_without_status, status_file -> [meta_without_status.id, status_file] },
            by: 0
        )
        .map { _id, meta, payload_file, data_files, status_file ->
            def updated_meta = meta.clone()
            
            // Read status from YAML file (simple text parsing)
            def status_text = status_file.text
            def process_status = status_text.contains('status: "SUCCESS"') ? 'SUCCESS' : 'FAILED'
            
            // Simple status: 'success' or 'failed'
            updated_meta.status = process_status == 'SUCCESS' ? 'success' : 'failed'
            
            return [ updated_meta, payload_file, data_files ]
        }
    
    // Step 2: Validate payloads (validation handles upstream failures internally)
    PAYLOAD_VALIDATE ( generation_with_status )
    ch_versions = ch_versions.mix(PAYLOAD_VALIDATE.out.versions.first())
    ch_status = ch_status.mix(PAYLOAD_VALIDATE.out.status)

    // Update meta.status after validation using status.yml
    validation_with_status = PAYLOAD_VALIDATE.out.validated_payload_files
        .map { meta, payload_file, data_files -> [meta.id, meta, payload_file, data_files] }
        .join(
            PAYLOAD_VALIDATE.out.status.map { meta_without_status, status_file -> [meta_without_status.id, status_file] },
            by: 0
        )
        .map { _id, meta, payload_file, data_files, status_file ->
            def updated_meta = meta.clone()
            
            // Read validation status from YAML file
            def status_text = status_file.text
            def validation_process_status = status_text.contains('status: "SUCCESS"') ? 'SUCCESS' : 'FAILED'
            
            // Only successful if validation passed AND previous status was success
            updated_meta.status = (validation_process_status == 'SUCCESS' && meta.status == 'success') ? 'success' : 'failed'
            
            return [ updated_meta, payload_file, data_files ]
        }

    // Simple filtering based on meta.status
    successful_analyses = validation_with_status
        .filter { meta, _payload_file, _data_files -> 
            meta.status == 'success' 
        }
    
    failed_analyses = validation_with_status
        .filter { meta, _payload_file, _data_files ->
            meta.status == 'failed'
        }

    emit:
    // Main outputs
    all_analyses     = validation_with_status       // channel: [ val(meta_with_status), payload.json, [data_files] ]
    successful       = successful_analyses          // channel: [ val(meta), payload.json, [data_files] ] where meta.status == 'success'
    failed           = failed_analyses              // channel: [ val(meta), payload.json, [data_files] ] where meta.status == 'failed'
    
    // Status and versions for full workflow aggregation
    status           = ch_status                    // channel: [ val(meta_without_status), status.yml ] - individual process status files
    versions         = ch_versions                  // channel: [ versions.yml ]
}
