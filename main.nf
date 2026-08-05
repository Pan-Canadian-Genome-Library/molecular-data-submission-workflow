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

    // Log active profiles, grouped by category
    def activeProfiles          = workflow.profile?.tokenize(',')*.trim() ?: []
    def knownEnvProfiles        = ['sd4h_dev', 'sd4h_staging', 'sd4h_qa', 'sd4h_prod', 'cumulus_dev', 'cumulus_qa']
    def knownStorageProfiles    = ['s3', 's3_compatible', 's3_anonymous']
    def knownContainerProfiles  = ['conda', 'mamba', 'docker', 'arm', 'singularity', 'apptainer', 'podman', 'shifter', 'charliecloud', 'wave']
    def knownSchedulerProfiles  = ['slurm']
    def knownOtherProfiles      = ['debug', 'test', 'test_sequential', 'test_parallel', 'test_full']

    def detectedEnv       = activeProfiles.findAll { it in knownEnvProfiles }
    def detectedStorage   = activeProfiles.any { it in knownStorageProfiles } ? ['s3'] : []
    def detectedContainer = activeProfiles.findAll { it in knownContainerProfiles }
    def detectedScheduler = activeProfiles.findAll { it in knownSchedulerProfiles }
    def detectedOther     = activeProfiles.findAll { it in knownOtherProfiles }

    // Scan the fileName column in file_metadata to detect all storage backends in use
    if (params.file_metadata) {
        try {
            def fmLines = file(params.file_metadata).readLines()
            if (fmLines.size() > 1) {
                def headers     = fmLines[0].split('\t')*.trim()
                def fileNameIdx = headers.findIndexOf { it == 'fileName' }
                if (fileNameIdx >= 0) {
                    def fileNames = fmLines.drop(1).collect { it.split('\t')[fileNameIdx] }
                    if (!detectedStorage.contains('s3')      && fileNames.any { it.startsWith('s3://') })                detectedStorage << 's3'
                    if (!detectedStorage.contains('azure')   && fileNames.any { it.startsWith('az://') })                detectedStorage << 'azure'
                    if (!detectedStorage.contains('gcs')     && fileNames.any { it.startsWith('gs://') })                detectedStorage << 'gcs'
                    if (!detectedStorage.contains('http')    && fileNames.any { it ==~ /^https?:\/\/.+/ })               detectedStorage << 'http'
                    if (!detectedStorage.contains('ftp')     && fileNames.any { it ==~ /^ftp:\/\/.+/ })                  detectedStorage << 'ftp'
                    if (!detectedStorage.contains('local')   && fileNames.any { !(it ==~ /^[a-z][a-z0-9+\-.]*:\/\//) }) detectedStorage << 'local'
                }
            }
        } catch (Exception e) {
            // file_metadata may not be locally readable at startup (e.g. itself on a remote store) — skip
        }
    }

    def allKnown     = knownEnvProfiles + knownStorageProfiles + knownContainerProfiles + knownSchedulerProfiles + knownOtherProfiles
    def unrecognized = activeProfiles.findAll { !(it in allKnown) }

    log.info "🚀 Active profiles:"
    if (detectedEnv) {
        log.info "   - Environment : ${detectedEnv.join(', ')}"
    } else {
        log.warn "   - Environment : ⚠️  No environment profile selected — use -profile <env> (e.g. sd4h_dev, sd4h_prod)"
    }
    log.info "   - File Storage: ${detectedStorage   ? detectedStorage.join(', ')   : 'Local'}"
    log.info "   - Container   : ${detectedContainer ? detectedContainer.join(', ') : 'None'}"
    log.info "   - Scheduler   : ${detectedScheduler ? detectedScheduler.join(', ') : 'Local'}"
    if (detectedOther)   log.info "   - Other       : ${detectedOther.join(', ')}"
    if (unrecognized)    log.warn "   ⚠️  Unrecognized profiles: ${unrecognized.join(', ')}"

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
    log.info "   - path_to_files_directory: ${params.path_to_files_directory ?: 'Not Provided'}"
    log.info "   - skip_upload: ${params.skip_upload}"
    log.info "   - allow_duplicates: ${params.allow_duplicates}"
    // Notify user if running in validation-only mode
    if (params.skip_upload) {
        log.info """
╔══════════════════════════════════════════════════════════════════════════════╗
║                     🧪 VALIDATION-ONLY MODE ENABLED                          ║
║         '--skip_upload' is enabled — no data will be submitted.              ║
║         The workflow will perform local validation steps only.               ║
║         To run a real submission, omit the '--skip_upload' flag              ║
║         and provide a valid token via '--token'.                             ║
╚══════════════════════════════════════════════════════════════════════════════╝
""".stripIndent()
    }
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
    if (params.token == null && params.skip_upload == false){
        startup_error_details.add("'token' was not provided, please provide the variable via the '--token' flag or in config.")
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
        params.path_to_files_directory ? channel.fromPath(params.path_to_files_directory) : []
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
