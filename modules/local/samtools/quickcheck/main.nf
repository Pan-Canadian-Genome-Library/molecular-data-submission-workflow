process SAMTOOLS_QUICKCHECK {
    tag "$meta.id:$alignment_file.baseName"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0':
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(payload), path(alignment_file), path(index_files)

    output:
    tuple val(meta), path(payload), path(alignment_file), path(index_files), emit: ch_validated_file
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
    QUICKCHECK_EXIT_CODE=0
    ERROR_DETAILS=""
    
    # Check if upstream process was successful
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed, skipping samtools quickcheck for ${alignment_file}"
        QUICKCHECK_EXIT_CODE=1
        ERROR_DETAILS="Skipped samtools quickcheck due to upstream failure"
    else
        echo "Running samtools quickcheck on: ${alignment_file}"
        samtools quickcheck -v "${alignment_file}"
        if [ \$? -ne 0 ]; then
            QUICKCHECK_EXIT_CODE=1
            # Capture error details from .command.err file
            if [ -f ".command.err" ] && [ -s ".command.err" ]; then
                # Remove ANSI color codes, carriage returns, and filter out empty lines
                ERROR_DETAILS=\$(sed 's/\\x1b\\[[0-9;]*m//g' ".command.err" | tr '\\r' '\\n' | grep -v '^[[:space:]]*\$' || cat ".command.err" | sed 's/\\x1b\\[[0-9;]*m//g')
            else
                ERROR_DETAILS="Alignment file ${alignment_file}: Failed samtools quickcheck"
            fi
            echo "ERROR: samtools quickcheck failed for ${alignment_file}"
        else
            echo "SUCCESS: samtools quickcheck passed for ${alignment_file}"
        fi
    fi
    
    # Create status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${alignment_file.baseName}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$QUICKCHECK_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$QUICKCHECK_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    details:
        analysis_id: "${meta.id}"
        tool: "samtools_quickcheck"
        file: "${alignment_file}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message if validation failed
    if [ \$QUICKCHECK_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "    error_details: |" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${alignment_file.baseName}_status.yml"
        echo "\$ERROR_DETAILS" | while IFS= read -r line; do
            echo "        \$line" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${alignment_file.baseName}_status.yml"
        done
    fi
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    
    # Exit handling
    if [ "${exit_on_error_str}" = "true" ] && [ \$QUICKCHECK_EXIT_CODE -ne 0 ]; then
        echo "samtools quickcheck failed and exit_on_error is true"
        exit \$QUICKCHECK_EXIT_CODE
    else
        echo "Continuing workflow (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"
    """
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_${alignment_file.baseName}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-08-11T10:30:00+00:00"
    details:
        analysis_id: "${meta.id}"
        tool: "samtools_quickcheck"
        file: "${alignment_file}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
