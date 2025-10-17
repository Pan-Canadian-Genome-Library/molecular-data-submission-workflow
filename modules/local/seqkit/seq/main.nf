process SEQKIT_SEQ {
    tag "$meta.id:$fastq_file.baseName"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    tuple val(meta), path(payload), path(fastq_file)

    output:
    tuple val(meta), path(payload), path(fastq_file), emit: ch_validated_file
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Initialize variables
    FASTQ_EXIT_CODE=0
    ERROR_DETAILS=""
    
    # Check if upstream process was successful
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed, skipping FASTQ validation for ${fastq_file}"
        FASTQ_EXIT_CODE=1
        ERROR_DETAILS="Skipped FASTQ validation due to upstream failure"
    else
        echo "Running seqkit seq validation on: ${fastq_file}"
        # Execute seqkit seq and capture stderr directly
        ERROR_OUTPUT=\$(seqkit seq --validate-seq --quiet "${fastq_file}" 2>&1)
        FASTQ_EXIT_CODE=\$?
        
        if [ \${FASTQ_EXIT_CODE} -ne 0 ]; then
            # Process the captured error output
            if [ -n "\${ERROR_OUTPUT}" ]; then
                # Remove ANSI color codes and clean up formatting
                ERROR_DETAILS=\$(echo "\${ERROR_OUTPUT}" | tr -d '\\033' | sed 's/\\[[0-9;]*m//g' | tr '\\r' '\\n' | grep -v '^[[:space:]]*\$' || echo "\${ERROR_OUTPUT}")
            else
                ERROR_DETAILS="FASTQ file ${fastq_file}: Failed seqkit seq validation"
            fi
            echo "ERROR: seqkit seq validation failed for ${fastq_file}"
        else
            echo "SUCCESS: seqkit seq validation passed for ${fastq_file}"
        fi
    fi
    
    # Create status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${fastq_file.baseName}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$FASTQ_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$FASTQ_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}"
        tool: "seqkit_seq"
        file: "${fastq_file}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message if validation failed
    if [ \$FASTQ_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "    error_details: |" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${fastq_file.baseName}_status.yml"
        echo "\$ERROR_DETAILS" | while IFS= read -r line; do
            echo "        \$line" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${fastq_file.baseName}_status.yml"
        done
    fi
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    
    # Exit handling
    if [ "${exit_on_error_str}" = "true" ] && [ \$FASTQ_EXIT_CODE -ne 0 ]; then
        echo "FASTQ validation failed and exit_on_error is true"
        exit \$FASTQ_EXIT_CODE
    else
        echo "Continuing workflow (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"
    """
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${fastq_file.baseName}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-08-11T10:30:00+00:00"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}"
        tool: "seqkit_seq"
        file: "${fastq_file}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.8.2
    END_VERSIONS
    """
}
