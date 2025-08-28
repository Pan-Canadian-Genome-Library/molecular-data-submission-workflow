process VALIDATION_METADATA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(payload), path(payload_files), path(specimen), path(sample), path(experiment), path(read_group)

    output:
    tuple val(meta), path(payload), path(payload_files), emit: ch_payload_files
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    def read_group_arg = (read_group && read_group != [] && read_group.toString() != '[]') ? "--read-group-file ${read_group}" : ""
    
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    METADATA_EXIT_CODE=0
    ERROR_DETAILS=""
    
    echo "Starting metadata validation..."
    echo "Analysis type: ${meta.type}"
    
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" = "failed" ]; then
        echo "Upstream process failed (meta.status=failed), skipping metadata validation"
        METADATA_EXIT_CODE=1
        ERROR_DETAILS="Skipped due to upstream failure"
    else
        echo "Upstream process successful or not set, proceeding with metadata validation"
        
        # Build Python script arguments - read_group is conditionally added based on availability
        PYTHON_ARGS="${payload}"
        PYTHON_ARGS="\$PYTHON_ARGS --specimen-file ${specimen}"
        PYTHON_ARGS="\$PYTHON_ARGS --sample-file ${sample}"
        PYTHON_ARGS="\$PYTHON_ARGS --experiment-file ${experiment}"
        PYTHON_ARGS="\$PYTHON_ARGS ${read_group_arg}"
        PYTHON_ARGS="\$PYTHON_ARGS --analysis-type ${meta.type}"
        
        # Run validation - Python script handles all file availability and requirement logic
        echo "Running validation with command: main.py \$PYTHON_ARGS"
        # Capture stderr for error details, let stdout flow normally (will be suppressed)
        main.py \$PYTHON_ARGS 2>validation_errors.tmp >/dev/null
        VALIDATION_EXIT_CODE=\$?
        
        if [ \$VALIDATION_EXIT_CODE -ne 0 ]; then
            METADATA_EXIT_CODE=1
            # Read error details from stderr (validation_errors.tmp)
            ERROR_DETAILS=\$(cat validation_errors.tmp)
        fi
        rm -f validation_errors.tmp
    fi
    
    # Create step-specific status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$METADATA_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$METADATA_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}"
        analysis_type: "${meta.type}"
        study_id: "${meta.study}"
        upstream_status: "${meta.status}"
        payload_file: "${payload}"
        read_group_file: "${read_group}"
        specimen_file: "${specimen}"
        sample_file: "${sample}"
        experiment_file: "${experiment}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message to status file if validation failed
    if [ \$METADATA_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        # Format multi-line error details properly for YAML
        echo "    error_details: |" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
        echo "\$ERROR_DETAILS" | sed 's/^/            /' >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    fi
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    
    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" = "true" ] && [ \$METADATA_EXIT_CODE -ne 0 ]; then
        echo "Metadata validation failed and exit_on_error is true, exiting with error"
        exit \$METADATA_EXIT_CODE
    else
        echo "Continuing workflow regardless of metadata validation result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"
    """
    # Create mock metadata status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-25T10:31:00+00:00"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}"
        analysis_type: "${meta.type}"
        study_id: "${meta.study}"
        upstream_status: "${meta.status}"
        payload_file: "${payload}"
        read_group_file: "${read_group}"
        specimen_file: "${specimen}"
        sample_file: "${sample}"
        experiment_file: "${experiment}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
