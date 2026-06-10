/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//include { paramsSummaryMap       } from 'plugin/nf-schema'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
//include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_molecular-data-submission-workflow_pipeline'

include {CHECK_SUBMISSION_DEPENDENCIES} from '../subworkflows/local/check_submission_dependencies/main.nf'
include {DATA_VALIDATION} from '../subworkflows/local/data_validation/main.nf'
include {DATA_UPLOAD} from '../subworkflows/local/data_upload/main.nf'
include {METADATA_PAYLOAD_GENERATION} from '../subworkflows/local/metadata_payload_generation/main.nf'
include {BATCH_RECEIPT_GENERATION} from '../subworkflows/local/batch_receipt_generation/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOLECULAR_DATA_SUBMISSION_WORKFLOW {

    take:
    file_metadata // Spreadsheet from --file_metadata
    analysis_metadata // Spreadsheet from --analysis_metadata
    workflow_metadata // Spreadsheet from --workflow_metadata
    read_group_metadata // Spreadsheet from --read_group_metadata 
    experiment_metadata // Spreadsheet from --experiment_metadata
    specimen_metadata // Spreadsheet from --specimen_metadata
    sample_metadata // Spreadsheet from --sample_metadata
    path_to_files_directory // file path to directory containing target files

    main:

    ch_versions = Channel.empty()
    ch_all_status = Channel.empty()  // Collect all [val(meta), *_status.yml] tuples from all processes
    
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/58
    CHECK_SUBMISSION_DEPENDENCIES(
        params.study_id, // string
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

    CHECK_SUBMISSION_DEPENDENCIES.out.status
    .flatten() // Flatten the collection into individual items
    .map { item ->
        // Extract meta and status_file from each item
        [item.meta, item.status_file]
    }
    .set { ch_status_tuples }

    ch_all_status = ch_all_status.mix(ch_status_tuples)

    // Debug: Check submission dependencies output channels
    if (params.debug_channels) {
        // log.info "🔍 CHECK_SUBMISSION_DEPENDENCIES outputs:"
        CHECK_SUBMISSION_DEPENDENCIES.out.submitted_analysis_channels.view { item ->
    "🔍 CHECK_SUBMISSION_DEPENDENCIES outputs: - submitted_analysis_channels - ID: ${item.meta.id}, status: ${item.meta.status}, analysis: ${item.analysis}, clinical: ${item.clinical}, files: ${item.files}"
}
    }

    // Manipulate the output channel from CHECK_SUBMISSION_DEPENDENCIES
    
    CHECK_SUBMISSION_DEPENDENCIES.out.submitted_analysis_channels
    .map { item ->
        def meta = item.meta
        def file_meta = item.analysis.files
        def analysis_meta = item.analysis.analysis
        def workflow_meta = item.analysis.workflow ?: []
        def data_files = item.files

        // Return the tuple in the required format
        [meta, file_meta, analysis_meta, workflow_meta, data_files]
    }.set{ ch_metadata_payload_input}

    CHECK_SUBMISSION_DEPENDENCIES.out.submitted_analysis_channels
    .map { item ->
        def meta = item.meta
        def specimen_meta = item.clinical.specimen ?: []
        def sample_meta = item.clinical.sample ?: []
        def experiment_meta = item.clinical.experiment ?: []
        def read_group_meta = item.clinical.read_group ?: []

        // Return the tuple in the required format
        [meta, [specimen_meta, sample_meta, experiment_meta, read_group_meta]]
    }.set{ ch_biospecimen_entity_input}

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/63

    METADATA_PAYLOAD_GENERATION(
        ch_metadata_payload_input
            .filter { meta, file_meta, analysis_meta, workflow_meta, data_files -> 
                (meta.status ?: 'pass') == 'pass' 
            } //channel: [val(meta), path(file_meta.tsv), path(analysis_meta.tsv), path(workflow_meta.tsv), [data_files]] per analysis
    )
    ch_versions = ch_versions.mix(METADATA_PAYLOAD_GENERATION.out.versions)
    ch_all_status = ch_all_status.mix(METADATA_PAYLOAD_GENERATION.out.status)

    // METADATA_PAYLOAD_GENERATION.out.all_analyses.subscribe{println "payload-generation: $it"}

    // Debug: Metadata payload generation output channels
    if (params.debug_channels) {
        // log.info "🔍 METADATA_PAYLOAD_GENERATION outputs:"
        METADATA_PAYLOAD_GENERATION.out.all_analyses.view { meta, payload, payload_files ->
            "🔍 METADATA_PAYLOAD_GENERATION outputs: - 📦 all_analyses - ID: ${meta.id}, status: ${meta.status}, payload: ${payload}, files: ${payload_files}"
        }
    }

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/60
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/61
    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/62

    // Join METADATA_PAYLOAD_GENERATION output with biospecimen entity data
    // Join by meta.id only, then combine status from both channels
    ch_joined_validation_input = METADATA_PAYLOAD_GENERATION.out.all_analyses
        .map { meta, payload, payload_files -> [meta.id, meta, payload, payload_files] } // Extract meta.id as join key
        .join(
            ch_biospecimen_entity_input
                .map { meta, entity_files -> [meta.id, meta, entity_files] }, // Extract meta.id as join key
            by: 0 // Join by meta.id (first element)
        )
        .map { id, meta_payload, payload, payload_files, meta_entity, entity_files ->
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
        .filter { meta, payload, payload_files, specimen, sample,experiment, read_group ->
            (meta.status ?: 'pass') == 'pass'
        }

    DATA_VALIDATION(
        ch_joined_validation_input // channel: [ val(meta), payload(json), [files], specimen(path), sample(path), experiment(path), read_group(path)]
    )
    ch_versions = ch_versions.mix(DATA_VALIDATION.out.versions)
    ch_all_status = ch_all_status.mix(DATA_VALIDATION.out.status)

    // Debug: Data validation output channels
    if (params.debug_channels) {
        // log.info "🔍 DATA_VALIDATION outputs:"
        DATA_VALIDATION.out.validated_payload_files.view { meta, payload, payload_files ->
            "🔍 DATA_VALIDATION outputs: - 📦 validated_payload_files - ID: ${meta.id}, status: ${meta.status}, payload: ${payload}, files: ${payload_files}"
        }
    }



    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/64
    ch_analysis = Channel.empty()
    if (!params.skip_upload) {
        DATA_UPLOAD(
            DATA_VALIDATION.out.validated_payload_files
                .filter { meta, _payload, _files -> meta.status == 'pass' }
        )
        ch_analysis = ch_analysis.mix(DATA_UPLOAD.out.analysis_json)
        ch_versions = ch_versions.mix(DATA_UPLOAD.out.versions)
        ch_all_status = ch_all_status.mix(DATA_UPLOAD.out.status)

        // Debug: Data upload output channels
        if (params.debug_channels) {
            // log.info "🔍 DATA_UPLOAD outputs:"
            DATA_UPLOAD.out.analysis_json.view { meta, analysis_json ->
                "🔍 DATA_UPLOAD outputs: - 📦 analysis_json - ID: ${meta.id}, status: ${meta.status}, analysis_json: ${analysis_json}"
            }
        }
    }   

    //https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/67
    // Generate batch receipts using all status files and analysis JSON
    // Group ch_all_status by meta after removing meta.status
    ch_grouped_status = ch_all_status
        .map { meta, status_file -> 
            // Remove meta.status before grouping to ensure proper grouping by meta
            def clean_meta = meta.clone()
            clean_meta.remove('status')
            clean_meta.remove('num_files')
            clean_meta.remove('num_validated_files')
            [clean_meta, status_file]
        }
        .groupTuple(by: 0) // Group by clean_meta (first element)
    
    // Left join the grouped status with analysis JSON - keep all grouped status records
    // When skip_upload is true, ch_analysis will be empty, so all records will have null analysis_json
    ch_status_analysis = ch_grouped_status
        .join(
            ch_analysis.map { meta, analysis_json -> 
                // Remove meta.status before joining
                def clean_meta = meta.clone()
                clean_meta.remove('status')
                clean_meta.remove('num_files')
                clean_meta.remove('num_validated_files')
                [clean_meta, analysis_json]
            },
            by: 0, // Join by clean_meta (first element)
            remainder: true // Left join - keep all records from ch_grouped_status
        )
        .map { meta, status_files, analysis_json -> 
            // analysis_json will be null for records that don't have matches in ch_analysis
            [meta, status_files, analysis_json ?: []]
        }

    BATCH_RECEIPT_GENERATION(ch_status_analysis)
    ch_versions = ch_versions.mix(BATCH_RECEIPT_GENERATION.out.versions)

    // Debug: Receipt generation output channels
    if (params.debug_channels) {
        BATCH_RECEIPT_GENERATION.out.batch_receipts.view { meta, batch_receipt_file_json, batch_receipt_file_tsv, total_analyses, successful_analyses, failed_analyses ->
            "🔍 BATCH_RECEIPT_GENERATION outputs: - 📦 batch_receipt_file - ID: ${meta.id}, batch_receipt_file_json: ${batch_receipt_file_json}, batch_receipt_file_tsv: ${batch_receipt_file_tsv}, total: ${total_analyses}, success: ${successful_analyses}, failed: ${failed_analyses}"
        }
    }

    //
    // Workflow completion notifications and messaging
    //
    BATCH_RECEIPT_GENERATION.out.batch_receipts
        .subscribe { meta, batch_receipt_file_json, batch_receipt_file_tsv, total_analyses, successful_analyses, failed_analyses ->

            // Print completion message with batch receipt location
            def receipt_path_json = batch_receipt_file_json.toAbsolutePath()
            def receipt_path_tsv = batch_receipt_file_tsv.toAbsolutePath()
            def total_count = total_analyses.text
            def success_count = successful_analyses.text
            def failed_count = failed_analyses.text
            def mode_label   = params.skip_upload ? "🧪 VALIDATION-ONLY MODE" : "🎉 WORKFLOW COMPLETED!"
            def success_label = params.skip_upload ? "✅ Passed validation:    " : "✅ Successful submissions:"
            def failed_label  = params.skip_upload ? "❌ Failed validation:    " : "❌ Failed submissions:    "
            log.info """
            ╔══════════════════════════════════════════════════════════════════════════════════╗
║                       ${mode_label}                      
            ╠══════════════════════════════════════════════════════════════════════════════════╣
            ║  Study: ${params.study_id}
            ║  Batch ID: ${meta.batch_id ?: meta.id}                                           
            ║  Total in this batch:    ${total_count}                                           
            ║  ${success_label} ${success_count}                                   
            ║  ${failed_label} ${failed_count}  
            ║                                      
            ║  📋 BATCH RECEIPT GENERATED:                                                     
            ║  📁 JSON RECEIPT Location: ${receipt_path_json}                                  
            ║  📁 TSV RECEIPT Location: ${receipt_path_tsv}                                   
            ║                                                                                  
            ║                                                                                  
            ║  ℹ️  The batch receipt contains:                                                 
            ║     • Summary of all processed analyses                                          
            ║     • Status of each ${params.skip_upload ? 'validation' : 'submission'} step                                             
            ║     • ${params.skip_upload ? 'Validation results for each analysis' : 'Upload results and analysis IDs'}                                            
            ║     • Error details for any failed analyses                                      
            ╚══════════════════════════════════════════════════════════════════════════════════╝
            """.stripIndent()
            
        }
        // .subscribe { meta, batch_receipt_file_json, batch_receipt_file_tsv, total_analyses, successful_analyses, failed_analyses ->

        //     // Print completion message with batch receipt location
        //     def receipt_path_json = batch_receipt_file_json.toAbsolutePath()
        //     def receipt_path_tsv = batch_receipt_file_tsv.toAbsolutePath()
        //     def total_count = total_analyses.text
        //     def success_count = successful_analyses.text
        //     def failed_count = failed_analyses.text
        //     log.info """
        //     ╔══════════════════════════════════════════════════════════════════════════════════╗
        //     ║                       🎉 WORKFLOW COMPLETED! 🎉                      
        //     ╠══════════════════════════════════════════════════════════════════════════════════╣
        //     ║  Study: ${params.study_id}
        //     ║  Batch ID: ${meta.batch_id ?: meta.id}                                           
        //     ║  Total in this batch:    ${total_count}                                           
        //     ║  ✅ Successful submissions: ${success_count}                                   
        //     ║  ❌ Failed submissions:     ${failed_count}  
        //     ║                                      
        //     ║  📋 BATCH RECEIPT GENERATED:                                                     
        //     ║  📁 JSON RECEIPT Location: ${receipt_path_json}                                  
        //     ║  📁 TSV RECEIPT Location: ${receipt_path_tsv}                                   
        //     ║                                                                                  
        //     ║                                                                                  
        //     ║  ℹ️  The batch receipt contains:                                                 
        //     ║     • Summary of all processed analyses                                          
        //     ║     • Status of each submission step                                             
        //     ║     • Upload results and analysis IDs                                            
        //     ║     • Error details for any failed analyses                                      
        //     ╚══════════════════════════════════════════════════════════════════════════════════╝
        //     """.stripIndent()
            
        // }

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
