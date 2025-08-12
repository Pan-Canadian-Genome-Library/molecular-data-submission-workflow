process SEQKIT_SEQ {
    tag "$meta.id:$fastq_file.baseName"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.2--h9ee0642_0':
        'biocontainers/seqkit:2.8.2--h9ee0642_0' }"

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
        # Use seqkit seq to validate FASTQ format
        seqkit seq --validate-seq --quiet "${fastq_file}" > /dev/null 2>&1
        if [ \$? -ne 0 ]; then
            FASTQ_EXIT_CODE=1
            ERROR_DETAILS="FASTQ file ${fastq_file}: Failed seqkit seq validation"
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
    details:
        analysis_id: "${meta.id}"
        tool: "seqkit_seq"
        file: "${fastq_file}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message if validation failed
    if [ \$FASTQ_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "    error_details: |" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${fastq_file.baseName}_status.yml"
        echo "\$ERROR_DETAILS" | sed 's/^/            /' >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${fastq_file.baseName}_status.yml"
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
