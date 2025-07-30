process PREPROCESS_SUBMISSION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(file_metadata), path(analysis_metadata), path(workflow_metadata), path(read_group_metadata), path(experiment_metadata), path(specimen_metadata), path(sample_metadata), path(files_directory)

    output:
    tuple val(meta), path("file_meta_${meta.id}.tsv"),     emit: file_meta
    tuple val(meta), path("analysis_meta_${meta.id}.tsv"), emit: analysis_meta
    tuple val(meta), path("workflow_meta_${meta.id}.tsv"), emit: workflow_meta
    tuple val(meta), path("specimen_${meta.id}.tsv"),     emit: specimen
    tuple val(meta), path("sample_${meta.id}.tsv"),       emit: sample
    tuple val(meta), path("experiment_${meta.id}.tsv"),   emit: experiment
    tuple val(meta), path("read_group_${meta.id}.tsv"),   emit: read_group
    tuple val(meta), path("data_files_${meta.id}.txt"),   emit: data_files
    tuple val(meta), path("*_status.yml"),                emit: status
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Initialize variables
    PREPROCESS_EXIT_CODE=0
    ERROR_DETAILS=""
    
    echo "Installing required Python packages..."
    
    # Create a temporary directory for package installation
    TEMP_PYTHON_LIB="\$(mktemp -d)/python_packages"
    mkdir -p "\$TEMP_PYTHON_LIB"
    
    # Install to temporary directory and add to Python path
    pip install --target "\$TEMP_PYTHON_LIB" PyYAML
    if [ \$? -ne 0 ]; then
        echo "Failed to install PyYAML package" >&2
        PREPROCESS_EXIT_CODE=1
        ERROR_DETAILS="Failed to install required Python dependencies"
    else
        # Set PYTHONPATH to include our temporary package directory
        export PYTHONPATH="\$TEMP_PYTHON_LIB:\${PYTHONPATH:-}"

        main.py \\
            "${meta.id}" \\
            "${meta.study}" \\
            "${meta.type}" \\
            "${meta.status}" \\
            "${file_metadata}" \\
            "${analysis_metadata}" \\
            "${workflow_metadata}" \\
            "${read_group_metadata}" \\
            "${experiment_metadata}" \\
            "${specimen_metadata}" \\
            "${sample_metadata}" \\
            "${files_directory}" 2>preprocess_errors.tmp

        PREPROCESS_EXIT_CODE=\$?
        
        # Capture error details if preprocessing failed
        if [ \$PREPROCESS_EXIT_CODE -ne 0 ] && [ -f "preprocess_errors.tmp" ] && [ -s "preprocess_errors.tmp" ]; then
            ERROR_DETAILS=\$(cat preprocess_errors.tmp)
        fi
        rm -f preprocess_errors.tmp
    fi
    
    # Create step-specific status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$PREPROCESS_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$PREPROCESS_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    details:
        analysis_id: "${meta.id}"
        analysis_type: "${meta.type}"
        study_id: "${meta.study}"
        upstream_status: "${meta.status}"
        file_metadata: "${file_metadata}"
        analysis_metadata: "${analysis_metadata}"
        workflow_metadata: "${workflow_metadata}"
        read_group_metadata: "${read_group_metadata}"
        experiment_metadata: "${experiment_metadata}"
        specimen_metadata: "${specimen_metadata}"
        sample_metadata: "${sample_metadata}"
        files_directory: "${files_directory}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message to status file if preprocessing failed
    if [ \$PREPROCESS_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        # Format multi-line error details properly for YAML
        echo "        error_details: |" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
        echo "\$ERROR_DETAILS" | sed 's/^/            /' >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    fi
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    
    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" = "true" ] && [ \$PREPROCESS_EXIT_CODE -ne 0 ]; then
        echo "Preprocessing failed and exit_on_error is true, exiting with error"
        exit \$PREPROCESS_EXIT_CODE
    else
        echo "Continuing workflow regardless of preprocessing result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"
    """
    touch file_meta_${meta.id}.tsv analysis_meta_${meta.id}.tsv workflow_meta_${meta.id}.tsv
    touch specimen_${meta.id}.tsv sample_${meta.id}.tsv experiment_${meta.id}.tsv read_group_${meta.id}.tsv
    touch data_files_${meta.id}.txt ${meta.id}_molecular_data_submission_workflow_check_submission_dependencies_preprocess_submission_status.yml
    touch versions.yml
    """
}
