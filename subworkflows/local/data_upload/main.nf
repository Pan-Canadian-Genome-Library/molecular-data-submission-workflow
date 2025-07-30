// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SONG_SUBMIT } from '../../../modules/local/song/submit/main'
include { SONG_MANIFEST } from '../../../modules/local/song/manifest/main'
include { SCORE_UPLOAD } from '../../../modules/local/score/upload/main'
include { SONG_PUBLISH } from '../../../modules/local/song/publish/main'
include { SONG_GETANALYSIS } from '../../../modules/local/song/getanalysis/main'

workflow DATA_UPLOAD {

    take:
    validated_payload_files // channel: [ val(meta), payload, files ]

    main:

    ch_versions = Channel.empty()

    // Upload files using local SONG/SCORE modules with status propagation
    // Each process will check meta.status and skip if upstream failed
    
    // 1. Create new analysis
    SONG_SUBMIT( validated_payload_files )
    ch_versions = ch_versions.mix(SONG_SUBMIT.out.versions)

    // 2. Generate file manifest for upload - update meta.status based on submit result
    ch_manifest_input = SONG_SUBMIT.out.analysis_files
        .join(SONG_SUBMIT.out.status, by: 0)
        .map { meta, analysis_file, payload, files, status_file ->
            // Read status from YAML file and update meta
            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()
            updated_meta.status = status_value
            [updated_meta, analysis_file, payload, files]
        }
    
    SONG_MANIFEST( ch_manifest_input )
    ch_versions = ch_versions.mix(SONG_MANIFEST.out.versions)

    // 3. Upload to SCORE - update meta.status based on manifest result  
    ch_upload_input = SONG_MANIFEST.out.manifest_upload
        .join(SONG_MANIFEST.out.status, by: 0)
        .map { meta, analysis_file, manifest, upload_files, status_file ->
            // Read status from YAML file and update meta
            def status_content = status_file.text  
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()
            updated_meta.status = status_value
            [updated_meta, analysis_file, manifest, upload_files]
        }
    
    SCORE_UPLOAD( ch_upload_input )
    ch_versions = ch_versions.mix(SCORE_UPLOAD.out.versions)

    // 4. Publish the analysis - update meta.status based on upload result
    ch_publish_input = SCORE_UPLOAD.out.ready_to_publish
        .join(SCORE_UPLOAD.out.status, by: 0)
        .map { meta, analysis_file, status_file ->
            // Read status from YAML file and update meta
            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass' 
            def updated_meta = meta.clone()
            updated_meta.status = status_value
            [updated_meta, analysis_file]
        }
    
    SONG_PUBLISH( ch_publish_input )
    ch_versions = ch_versions.mix(SONG_PUBLISH.out.versions)

    // 5. Get analysis details from server - update meta.status based on publish result
    ch_getanalysis_input = SONG_PUBLISH.out.analysis_id
        .join(SONG_PUBLISH.out.status, by: 0)
        .map { meta, analysis_file, status_file ->
            // Read status from YAML file and update meta
            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()
            updated_meta.status = status_value
            [updated_meta, analysis_file]
        }
    
    SONG_GETANALYSIS( ch_getanalysis_input )
    ch_versions = ch_versions.mix(SONG_GETANALYSIS.out.versions)

    // Collect all status files from upload pipeline
    ch_all_status = Channel.empty()
    ch_all_status = ch_all_status.mix(SONG_SUBMIT.out.status)
    ch_all_status = ch_all_status.mix(SONG_MANIFEST.out.status)
    ch_all_status = ch_all_status.mix(SCORE_UPLOAD.out.status)
    ch_all_status = ch_all_status.mix(SONG_PUBLISH.out.status)
    ch_all_status = ch_all_status.mix(SONG_GETANALYSIS.out.status)

    emit:
    analysis_json = SONG_GETANALYSIS.out.analysis_json     // channel: [ val(meta), path(*.analysis.json) ]
    status        = ch_all_status                           // channel: [ val(meta), path(*_status.yml) ]
    versions      = ch_versions                             // channel: [ versions.yml ]
}
