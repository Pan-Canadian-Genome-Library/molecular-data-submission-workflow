/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//include { paramsSummaryMap       } from 'plugin/nf-schema'

//include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
//include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_molecular-data-submission-workflow_pipeline'

include {CHECK_SUBMISSION_DEPENDENCIES} from '../subworkflows/local/check_submission_dependencies.nf'
include {CLINICAL_SERVICE_DATA_SUBMISSION} from '../subworkflows/local/clinical_service_data_submission.nf'
include {DATA_VALIDATION} from '../subworkflows/local/data_validation.nf'
include {METADATA_PAYLOAD_GENERATION} from '../subworkflows/local/metadata_payload_generation.nf'
include {SKIP_DUPLICATE_SONG_SCORE_UPLOAD} from '../subworkflows/local/skip_duplicate_song_score_upload.nf'
include {SUBMISSION_RECEIPT} from '../subworkflows/local/submission_receipt.nf'
include {SONG_SCORE_UPLOAD} from '../subworkflows/icgc-argo-workflows/song_score_upload/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOLECULAR_DATA_SUBMISSION_WORKFLOW {

    take:
    study_id
    token
    file_metadata // Spreadsheet from --file_metadata
    analysis_metadata // Spreadsheet from --analysis_metadata
    workflow_metadata // Spreadsheet from --workflow_metadata
    read_group_metadata // Spreadsheet from --read_group_metadata 
    experiment_metadata // Spreadsheet from --experiment_metadata
    specimen_metadata // Spreadsheet from --specimen_metadata
    sample_metadata // Spreadsheet from --sample_metadata
    path_to_files_directory // file path to directory containing target files
    skip_duplicate_check // pipeline flag from --skip_duplicate_check
    skip_upload // pipeline flag from --skip_duplicate_check
    main:

    ch_versions = Channel.empty()
    
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/58
    CHECK_SUBMISSION_DEPENDENCIES(
        study_id,
        token,
        file_metadata, // Spreadsheet
        analysis_metadata, // Spreadsheet
        workflow_metadata, // Spreadsheet
        read_group_metadata, // Spreadsheet
        experiment_metadata, // Spreadsheet
        specimen_metadata, // Spreadsheet
        sample_metadata, // Spreadsheet
        path_to_files_directory // file path
    )
    ch_versions = ch_versions.mix(CHECK_SUBMISSION_DEPENDENCIES.out.versions)

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/59
    CLINICAL_SERVICE_DATA_SUBMISSION(
        CHECK_SUBMISSION_DEPENDENCIES.out.clinical_upload // [ val(meta), [csv] ] 
    )

    // In case no health data needs to be submitted
    ch_unsuccessful_registeration = Channel.empty()
    ch_successful_registeration = Channel.empty()   
    ch_unsuccessful_registeration = ch_unsuccessful_registeration.concat(CLINICAL_SERVICE_DATA_SUBMISSION.out.unsuccessful_registeration)
    ch_successful_registeration = ch_successful_registeration.concat(CLINICAL_SERVICE_DATA_SUBMISSION.out.successful_registeration)

    ch_versions = ch_versions.mix(CLINICAL_SERVICE_DATA_SUBMISSION.out.versions)

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/63
    METADATA_PAYLOAD_GENERATION(
        CHECK_SUBMISSION_DEPENDENCIES.out.molecular_files_to_upload // channel [ val(meta) [files] ] per analysis
    )
    ch_versions = ch_versions.mix(METADATA_PAYLOAD_GENERATION.out.versions)

    METADATA_PAYLOAD_GENERATION.out.files_to_upload.subscribe{println "TEST1 $it"}

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/60
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/61
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/62

    DATA_VALIDATION(
        METADATA_PAYLOAD_GENERATION.out.files_to_upload // channel: [ val(meta), JSON payload, [genomic files] ] per payload 
    )
    ch_versions = ch_versions.mix(DATA_VALIDATION.out.versions)

    DATA_VALIDATION.out.successful_validation.subscribe{println "TEST2A $it"}

    DATA_VALIDATION.out.unsuccessful_validation.subscribe{println "TEST2B $it"}

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/64
    if (!skip_upload){ 
        if (skip_duplicate_check){
            // SKIP_DUPLICATE_SONG_SCORE_UPLOAD(
            DATA_VALIDATION.out.successful_validation
            unsuccessful_upload=[]
            successful_upload=[]
            // )
            //ch_versions = ch_versions.mix(SKIP_DUPLICATE_SONG_SCORE_UPLOAD.out.versions)
        } else {
            unsuccessful_upload=[]
            successful_upload=[]
            //SONG_SCORE_UPLOAD(DATA_VALIDATION.out.successful_validation.map{ meta, payload, files , checkA , checkB , checkC -> [meta , payload, files] })
            //ch_versions = ch_versions.mix(SONG_SCORE_UPLOAD.out.versions)
        }
    }

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/67
    SUBMISSION_RECEIPT(
       CHECK_SUBMISSION_DEPENDENCIES.out.unsuccessful_dependency, // [val(meta),[csv]]
       ch_unsuccessful_registeration, // channel: [ val(meta), [csv] ] multiple CSVs per entity
       DATA_VALIDATION.out.unsuccessful_validation,  //[ val(meta), JSON payload, [genomic files] , metadata sanity boolean, sequencing integrity boolean , cross validation boolean]
       unsuccessful_upload,
       CLINICAL_SERVICE_DATA_SUBMISSION.out.successful_registeration, // channel: [ val(meta), [csv] ] multiple CSVs per entity
       ch_successful_registeration,  //[ val(meta), JSON payload, [genomic files] , metadata sanity boolean, sequencing integrity boolean , cross validation boolean]
       successful_upload
    )
    ch_versions = ch_versions.mix(SUBMISSION_RECEIPT.out.versions)

    //
    // Collate and save software versions
    //
    // softwareVersionsToYAML(ch_versions)
    //     .collectFile(
    //         storeDir: "${params.outdir}/pipeline_info",
    //         name:  ''  + 'pipeline_software_' +  ''  + 'versions.yml',
    //         sort: true,
    //         newLine: true
    //     ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
