#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pan-Canadian-Genome-Library/molecular-data-submission-workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MOLECULAR_DATA_SUBMISSION_WORKFLOW  } from './workflows/molecular-data-submission-workflow'
//include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_molecular-data-submission-workflow_pipeline'
//include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_molecular-data-submission-workflow_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow PANCANADIANGENOMELIBRARY_MOLECULAR_DATA_SUBMISSION_WORKFLOW {

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
    skip_upload // pipeline flag from --skip_upload

    main:

    //
    // WORKFLOW: Run pipeline
    //
    MOLECULAR_DATA_SUBMISSION_WORKFLOW (
        study_id,
        token,
        file_metadata,
        analysis_metadata,
        workflow_metadata,
        read_group_metadata,
        experiment_metadata,
        specimen_metadata,
        sample_metadata,
        path_to_files_directory,
        skip_duplicate_check,
        skip_upload
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    // PIPELINE_INITIALISATION (
    //     params.version,
    //     params.validate_params,
    //     params.monochrome_logs,
    //     args,
    //     params.outdir,
    //     params.input
    // )

    //
    // WORKFLOW: Run main workflow
    //
    PANCANADIANGENOMELIBRARY_MOLECULAR_DATA_SUBMISSION_WORKFLOW (
        params.study_id,
        params.token,
        params.file_metadata,
        params.analysis_metadata,
        params.workflow_metadata,
        params.read_group_metadata,
        params.experiment_metadata,
        params.specimen_metadata,
        params.sample_metadata,
        params.path_to_files_directory,
        params.skip_duplicate_check,
        params.skip_upload
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    // PIPELINE_COMPLETION (
    //     params.outdir,
    //     params.monochrome_logs,
        
        
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
