

process CLINICAL_SUBMISSION {
    tag "$meta.id"
    label 'process_single'
    maxForks 1

container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:2.0.1--pyhdfd78af_0':
        'biocontainers/ampcombi:2.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(analysis), val(clinical), val(files), path(status_file), path(relational_mapping), path(analysis_types), path(data_directory)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), val(analysis), val(clinical), val(files) , path("*clinical_submission_status.yml") , emit : status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error.toString()
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_file = clinical.sample!=null ? "--sample_metadata ${clinical.sample}" : ""
    def specimen_file = clinical.specimen!=null ? "--specimen_metadata ${clinical.specimen}" : ""
    def experiment_file = clinical.experiment!=null ? "--experiment_metadata ${clinical.experiment}" : ""
    def read_group_file = clinical.read_group!=null ? "--read_group_metadata ${clinical.read_group}" : ""
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass' }" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping payload generation"
        GENERATION_EXIT_CODE=1
        ERROR_DETAILS="Submit clinical generation due to upstream failure"

    else
        echo "Upstream process successful, proceeding with payload generation"
        
        # Run main.py once and capture both exit code and error output
        main.py \\
        ${sample_file} \
        ${specimen_file} \
        ${experiment_file} \
        ${read_group_file} \
        --clinical_url ${params.clinical_url} \
        --relational_mapping ${relational_mapping} \
        --token ${params.token} \
        --study_id ${meta.study} 2>generation_errors.tmp

        GENERATION_EXIT_CODE=\$?

        # Capture error details if generation failed
        ERROR_DETAILS=""
        if [ \$GENERATION_EXIT_CODE -ne 0 ]; then
            # Read all error details from captured stderr
            if [ -f "generation_errors.tmp" ] && [ -s "generation_errors.tmp" ]; then
                ERROR_DETAILS=\$(cat "generation_errors.tmp")
            else
                ERROR_DETAILS="Script execution failed"
            fi
        fi
    fi

    # Create step-specific status file (meta.id used for grouping in channel)
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$GENERATION_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$GENERATION_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}" 
        sample_meta : "${clinical.sample}"
        specimen_meta :  "${clinical.specimen}"
        experiment_meta : "${clinical.experiment}"
        read_group_meta : "${clinical.read_group}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message to status file if generation failed
    if [ \$GENERATION_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "    error_details: |" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
        echo "\$ERROR_DETAILS" | while IFS= read -r line; do
            echo "        \$line" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
        done
    fi
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    
    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" == "true" ] && [ \$GENERATION_EXIT_CODE -ne 0 ]; then
        echo "Payload generation failed and exit_on_error is true, exiting with error"
        exit \$GENERATION_EXIT_CODE
    else
        echo "Continuing workflow regardless of generation result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error.toString()
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_file = clinical.sample!=null ? "--sample_metadata ${clinical.sample}" : ""
    def specimen_file = clinical.specimen!=null ? "--specimen_metadata ${clinical.specimen}" : ""
    def experiment_file = clinical.experiment!=null ? "--experiment_metadata ${clinical.experiment}" : ""
    def read_group_file = clinical.read_group!=null ? "--read_group_metadata ${clinical.read_group}" : ""
    """
    # Create mock step-specific status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-22T10:30:00+00:00"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}" 
        sample_meta : "${clinical.sample}"
        specimen_meta :  "${clinical.specimen}"
        experiment_meta : "${clinical.experiment}"
        read_group_meta : "${clinical.read_group}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}