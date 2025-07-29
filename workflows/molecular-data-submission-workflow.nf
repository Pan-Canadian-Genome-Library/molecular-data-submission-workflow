/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//include { paramsSummaryMap       } from 'plugin/nf-schema'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
//include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_molecular-data-submission-workflow_pipeline'

include {CHECK_SUBMISSION_DEPENDENCIES} from '../subworkflows/local/check_submission_dependencies/main.nf'
// include {CLINICAL_SERVICE_DATA_SUBMISSION} from '../subworkflows/local/clinical_service_data_submission.nf'
include {DATA_VALIDATION} from '../subworkflows/local/data_validation/main.nf'
include {DATA_UPLOAD} from '../subworkflows/local/data_upload/main.nf'
include {METADATA_PAYLOAD_GENERATION} from '../subworkflows/local/metadata_payload_generation/main.nf'
include {BATCH_RECEIPT_GENERATION} from '../subworkflows/local/batch_receipt_generation/main.nf'
//include {SKIP_DUPLICATE_SONG_SCORE_UPLOAD} from '../subworkflows/local/skip_duplicate_song_score_upload.nf'

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
    ch_all_status = Channel.empty()  // Collect all [val(meta), *_status.yml] tuples from all processes
    
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/58
    CHECK_SUBMISSION_DEPENDENCIES(
        file_metadata, // Spreadsheet
        analysis_metadata, // Spreadsheet
        workflow_metadata, // Spreadsheet
        read_group_metadata, // Spreadsheet
        experiment_metadata, // Spreadsheet
        specimen_metadata, // Spreadsheet
        sample_metadata, // Spreadsheet
        path_to_files_directory // file path to directory containing target files
    )
    ch_versions = ch_versions.mix(CHECK_SUBMISSION_DEPENDENCIES.out.versions)
    ch_all_status = ch_all_status.mix(CHECK_SUBMISSION_DEPENDENCIES.out.status)

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/63
    METADATA_PAYLOAD_GENERATION(
        CHECK_SUBMISSION_DEPENDENCIES.out.molecular_files_to_upload
            .filter { meta, _file_meta, _analysis_meta, _workflow_meta, _data_files -> 
                (meta.status ?: 'pass') == 'pass' 
            } //channel: [val(meta), path(file_meta.tsv), path(analysis_meta.tsv), path(workflow_meta.tsv), [data_files]] per analysis
    )
    ch_versions = ch_versions.mix(METADATA_PAYLOAD_GENERATION.out.versions)
    ch_all_status = ch_all_status.mix(METADATA_PAYLOAD_GENERATION.out.status)

    METADATA_PAYLOAD_GENERATION.out.all_analyses.subscribe{println "payload-generation: $it"}

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/60
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/61
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/62

    // Join METADATA_PAYLOAD_GENERATION output with biospecimen entity data
    // Join by meta.id only, then combine status from both channels
    ch_joined_validation_input = METADATA_PAYLOAD_GENERATION.out.all_analyses
        .map { meta, payload, payload_files -> [meta.id, meta, payload, payload_files] } // Extract meta.id as join key
        .join(
            CHECK_SUBMISSION_DEPENDENCIES.out.biospecimen_entity
                .map { meta, entity_files -> [meta.id, meta, entity_files] }, // Extract meta.id as join key
            by: 0 // Join by meta.id (first element)
        )
        .map { _id, meta_payload, payload, payload_files, meta_entity, entity_files ->
            // Now we have both meta objects separately
            // Check status from both channels
            def payload_status = meta_payload.status ?: 'pass'
            def entity_status = meta_entity.status ?: 'pass'
            
            // If either upstream process failed, mark as failed
            def combined_status = (payload_status == 'failed' || entity_status == 'failed') ? 'failed' : 'pass'
            
            // Update meta with combined status, keeping payload meta as base
            def updated_meta = meta_payload.clone()
            updated_meta.status = combined_status
            
            // Extract individual entity files from the entity_files collection
            // Assuming entity_files contains [specimen, sample, experiment, read_group] files
            def specimen = entity_files[0]
            def sample = entity_files[1] 
            def experiment = entity_files[2]
            def read_group = entity_files[3]
            
            return [updated_meta, payload, payload_files, specimen, sample, experiment, read_group]
        }
        .filter { meta, _payload, _payload_files, _specimen, _sample, _experiment, _read_group ->
            (meta.status ?: 'pass') == 'pass'
        }

    DATA_VALIDATION(
        ch_joined_validation_input // channel: [ val(meta), payload(json), [files], specimen(path), sample(path), experiment(path), read_group(path)]
    )
    ch_versions = ch_versions.mix(DATA_VALIDATION.out.versions)
    ch_all_status = ch_all_status.mix(DATA_VALIDATION.out.status)

    DATA_VALIDATION.out.validated_payload_files.subscribe{println "data-validated: $it"}


    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/64
    ch_analysis = Channel.empty()
    if (!skip_upload) {
        DATA_UPLOAD(
            DATA_VALIDATION.out.validated_payload_files
                .filter { meta, _payload, _files -> meta.status == 'pass' }
        )
        ch_analysis = ch_analysis.mix(DATA_UPLOAD.out.analysis_json)
        ch_versions = ch_versions.mix(DATA_UPLOAD.out.versions)
        ch_all_status = ch_all_status.mix(DATA_UPLOAD.out.status)
    }   

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/67
    // Generate batch receipts using all status files and analysis JSON
    // Group ch_all_status by meta after removing meta.status
    ch_grouped_status = ch_all_status
        .map { meta, status_file -> 
            // Remove meta.status before grouping to ensure proper grouping by meta
            def clean_meta = meta.clone()
            clean_meta.remove('status')
            [clean_meta, status_file]
        }
        .groupTuple(by: 0) // Group by clean_meta (first element)
        
    ch_grouped_status.subscribe { println "Grouped status: $it" }
    
    // Left join the grouped status with analysis JSON - keep all grouped status records
    // When skip_upload is true, ch_analysis will be empty, so all records will have null analysis_json
    ch_status_analysis = ch_grouped_status
        .join(
            ch_analysis.map { meta, analysis_json -> 
                // Remove meta.status before joining
                def clean_meta = meta.clone()
                clean_meta.remove('status')
                [clean_meta, analysis_json]
            },
            by: 0, // Join by clean_meta (first element)
            remainder: true // Left join - keep all records from ch_grouped_status
        )
        .map { meta, status_files, analysis_json -> 
            // analysis_json will be null for records that don't have matches in ch_analysis
            [meta, status_files, analysis_json ?: []]
        }
    ch_status_analysis.subscribe{println "status_analysis: $it"}
    BATCH_RECEIPT_GENERATION(ch_status_analysis)
    ch_versions = ch_versions.mix(BATCH_RECEIPT_GENERATION.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  ''  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_collated_versions                 // channel: [ path(versions.yml) ]
    all_status     = ch_all_status              // channel: [ val(meta), path(*_status.yml) ] - all status files from all processes

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
