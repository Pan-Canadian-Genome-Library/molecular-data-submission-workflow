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
workflow PCGL {

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

    //
    // WORKFLOW: Run pipeline
    //
    MOLECULAR_DATA_SUBMISSION_WORKFLOW (
        file_metadata,
        analysis_metadata,
        workflow_metadata,
        read_group_metadata,
        experiment_metadata,
        specimen_metadata,
        sample_metadata,
        path_to_files_directory
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

    // Log pipeline parameters
    log.info "🔧 Input parameters:"
    log.info "   - study_id: ${params.study_id}"
    log.info "   - analysis_metadata: ${params.analysis_metadata}"
    log.info "   - file_metadata: ${params.file_metadata}"
    log.info "   - workflow_metadata: ${params.workflow_metadata ?: 'Not Provided'}"
    log.info "   - read_group_metadata: ${params.read_group_metadata ?: 'Not Provided'}"
    log.info "   - experiment_metadata: ${params.experiment_metadata ?: 'Not Provided'}"
    log.info "   - specimen_metadata: ${params.specimen_metadata ?: 'Not Provided'}"
    log.info "   - sample_metadata: ${params.sample_metadata ?: 'Not Provided'}"
    log.info "   - path_to_files_directory: ${params.path_to_files_directory}"
    log.info "   - skip_upload: ${params.skip_upload}"
    log.info "   - allow_duplicates: ${params.allow_duplicates}"

    // Catch if dependencies are missing:
    startup_error_details=[]

    if (params.study_id == null){
         startup_error_details.add("'study_id' was not provided, please provide the variable via the '--study_id' flag or in config.")
    }
    if (params.analysis_metadata == null){
         startup_error_details.add("'analysis_metadata' was not provided, please provide the variable via the '--analysis_metdata' flag or in config.")
    }
    if (params.file_metadata == null){
         startup_error_details.add("'file_metadata' was not provided, please provide the variable via the '--file_metadata' flag or in config.")
    }
    if (params.token == null){
         startup_error_details.add("'token' was not provided, please provide the variable via the '--token' flag or in config.")
    }
    if (params.path_to_files_directory == null){
         startup_error_details.add("'path_to_files_directory' was not provided, please provide the variable via the '--path_to_files_directory' flag or in config.")
    }

    if (startup_error_details.size()!=0){
System.err.println """
╔══════════════════════════════════════════════════════════════════════════════╗
║                        🚨 WORKFLOW STOPPED                                   ║
║              Minimum requirements for data submission not met!               ║
╚══════════════════════════════════════════════════════════════════════════════╝
""".stripIndent()
System.err.println """
🔍 Issues found:
--------------------------------------------------------------------------------
""".stripIndent()
System.err.println """
${startup_error_details.join("\n")}
""".stripIndent()
System.err.println """
--------------------------------------------------------------------------------

""".stripIndent()
System.err.println """
Please fix the above issues and re-run the workflow.
                    """.stripIndent()
                    
                    // Flush the error stream before exiting
                    System.err.flush()

                    // Add a small delay to ensure output is displayed
                    Thread.sleep(100)

                    // Gracefully exit the workflow
                    System.exit(1)
    }
    //
    // WORKFLOW: Run main workflow
    //
    PCGL (
        channel.fromPath(params.file_metadata),
        channel.fromPath(params.analysis_metadata),
        params.workflow_metadata ? channel.fromPath(params.workflow_metadata) : [],
        params.read_group_metadata ? channel.fromPath(params.read_group_metadata) : [],
        params.experiment_metadata ? channel.fromPath(params.experiment_metadata) : [],
        params.specimen_metadata ? channel.fromPath(params.specimen_metadata) : [],
        params.sample_metadata ? channel.fromPath(params.sample_metadata) : [],
        channel.fromPath(params.path_to_files_directory)
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
