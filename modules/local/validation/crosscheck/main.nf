process VALIDATION_CROSSCHECK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(payload), path(payload_files)

    output:
    tuple val(meta), path(payload), path(payload_files), emit: ch_payload_files
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Initialize variables
    CROSSCHECK_EXIT_CODE=0
    ERROR_DETAILS=""
    
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping MD5 checksum validation"
        CROSSCHECK_EXIT_CODE=1
        ERROR_DETAILS="Skipped MD5 checksum validation due to upstream failure"
    elif grep -q '"error".*"payload_generation_failed"' "${payload}" 2>/dev/null; then
        echo "Detected placeholder payload file, skipping MD5 checksum validation"
        CROSSCHECK_EXIT_CODE=1
        ERROR_DETAILS="Skipped MD5 checksum validation - placeholder payload file detected"
    else
        echo "Upstream process successful, proceeding with MD5 checksum validation"
        
        # Run MD5 validation using external Python script
        # Capture stderr for error details, let stdout flow normally
        # Note: payload_files may contain multiple files, so we don't quote it to allow expansion
        main.py "${payload}" ${payload_files} 2>crosscheck_errors.tmp >/dev/null
        VALIDATION_EXIT_CODE=\$?
        
        if [ \$VALIDATION_EXIT_CODE -ne 0 ]; then
            CROSSCHECK_EXIT_CODE=1
            # Read error details from stderr
            ERROR_DETAILS="\${ERROR_DETAILS}\$(cat crosscheck_errors.tmp)"
        fi
        rm -f crosscheck_errors.tmp
    fi
    
    # Create step-specific status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$CROSSCHECK_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$CROSSCHECK_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}"
        payload_file: "${payload}"
        files_validated: "${payload_files}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message to status file if validation failed
    if [ \$CROSSCHECK_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
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
    if [ "${exit_on_error_str}" = "true" ] && [ \$CROSSCHECK_EXIT_CODE -ne 0 ]; then
        echo "MD5 checksum validation failed and exit_on_error is true, exiting with error"
        exit \$CROSSCHECK_EXIT_CODE
    else
        echo "Continuing workflow regardless of MD5 checksum validation result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"
    """
    # Create mock cross-check status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-25T10:32:00+00:00"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}"
        payload_file: "${payload}"
        files_validated: "${payload_files}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
