// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SONG_SUBMIT } from '../../../modules/local/song/submit/main'
include { SONG_MANIFEST } from '../../../modules/local/song/manifest/main'
include { SCORE_UPLOAD } from '../../../modules/local/score/upload/main'
include { SONG_PUBLISH } from '../../../modules/local/song/publish/main'

workflow DATA_UPLOAD {

    take:
    validated_payload_files // channel: [ val(meta), payload, files ]

    main:

    ch_versions = Channel.empty()

    // Upload files using local SONG/SCORE modules
    // Each process will handle its own error conditions and create appropriate status files
    
    // 1. Create new analysis
    SONG_SUBMIT( validated_payload_files )
    ch_versions = ch_versions.mix(SONG_SUBMIT.out.versions)

    // 2. Generate file manifest for upload
    SONG_MANIFEST( SONG_SUBMIT.out.analysis_files )
    ch_versions = ch_versions.mix(SONG_MANIFEST.out.versions)

    // 3. Upload to SCORE
    SCORE_UPLOAD( SONG_MANIFEST.out.manifest_upload )
    ch_versions = ch_versions.mix(SCORE_UPLOAD.out.versions)

    // 4. Publish the analysis
    SONG_PUBLISH( SCORE_UPLOAD.out.ready_to_publish )
    ch_versions = ch_versions.mix(SONG_PUBLISH.out.versions)

    // Collect all status files from upload pipeline
    ch_all_status = Channel.empty()
    ch_all_status = ch_all_status.mix(SONG_SUBMIT.out.status)
    ch_all_status = ch_all_status.mix(SONG_MANIFEST.out.status)
    ch_all_status = ch_all_status.mix(SCORE_UPLOAD.out.status)
    ch_all_status = ch_all_status.mix(SONG_PUBLISH.out.status)

    emit:
    analysis_id  = SONG_SUBMIT.out.analysis_id           // channel: [ val(meta), path(analysis_id.txt) ]
    status       = ch_all_status                        // channel: [ val(meta), path(*_status.yml) ]
    versions     = ch_versions                          // channel: [ versions.yml ]
}
