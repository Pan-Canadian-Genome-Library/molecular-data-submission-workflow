
process SCORE_UPLOAD {
    tag "$meta.id"
    label 'process_medium'

    container "${ params.file_transfer_container ?: 'ghcr.io/pan-canadian-genome-library/file-transfer' }:${ params.file_transfer_container_tag }"
    containerOptions "-v \$(pwd):/score-client/logs"

    input:
    tuple val(meta), path(analysis_id_file), path(manifest), path(upload)

    output:
    tuple val(meta), path("out/analysis_id.txt"),    emit: ready_to_publish
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    def file_manager_url = params.file_manager_url_upload ?: params.file_manager_url
    def file_transfer_url = params.file_transfer_url_upload ?: params.file_transfer_url
    def transport_parallel = params.file_transfer_transport_parallel ?: task.cpus
    def transport_mem = params.file_transfer_transport_mem ?: "2"
    def accessToken = params.token
    def VERSION = params.file_transfer_container_tag ?: 'edge'
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"

    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Create output directory first
    mkdir -p out

    # Copy analysis_id file to output for downstream processes  
    cp ${analysis_id_file} ./out/analysis_id.txt
    
    # Initialize variables
    UPLOAD_EXIT_CODE=0
    ERROR_DETAILS=""
    ANALYSIS_ID=""
    
    # Read analysis ID from file first (needed for both success and failure cases)
    if [ -f "${analysis_id_file}" ]; then
        ANALYSIS_ID=\$(cat ${analysis_id_file} | tr -d '\\n')
    fi
    
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping SCORE upload"
        UPLOAD_EXIT_CODE=1
        ERROR_DETAILS="Skipped SCORE upload due to upstream failure"
    else
        echo "Upstream process successful, proceeding with SCORE upload"
        
        export ACCESSTOKEN=${accessToken}
        export METADATA_URL=${file_manager_url}
        export STORAGE_URL=${file_transfer_url}
        export TRANSPORT_PARALLEL=${transport_parallel}
        export TRANSPORT_MEM=${transport_mem}

        # Execute SCORE upload
        score-client upload --manifest ${manifest} $args
        UPLOAD_EXIT_CODE=\${?}
        
        if [ \${UPLOAD_EXIT_CODE} -ne 0 ]; then
            # Capture error details, converting carriage returns to newlines and filtering for ERROR lines
            if [ -f ".command.err" ] && [ -s ".command.err" ]; then
                # Convert carriage returns to newlines, then get lines containing ERROR
                ERROR_DETAILS=\$(tr '\\r' '\\n' < ".command.err" | grep "ERROR" || echo "No ERROR lines found")
            else
                ERROR_DETAILS="Script execution failed - no error details available"
            fi
        fi
    fi

    # Create step-specific status file
    if [ \${UPLOAD_EXIT_CODE} -eq 0 ]; then
        STATUS_RESULT="SUCCESS"
    else
        STATUS_RESULT="FAILED"
    fi
    
    echo "process: \\"${task.process}\\"" > "${status_file_name}"
    echo "status: \\"\${STATUS_RESULT}\\"" >> "${status_file_name}"
    echo "exit_code: \${UPLOAD_EXIT_CODE}" >> "${status_file_name}"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${status_file_name}"
    echo "details:" >> "${status_file_name}"
    echo "    analysis_id: \\"\${ANALYSIS_ID}\\"" >> "${status_file_name}"
    echo "    file_manager_url: \\"${file_manager_url}\\"" >> "${status_file_name}"
    echo "    file_transfer_url: \\"${file_transfer_url}\\"" >> "${status_file_name}"
    echo "    manifest_file: \\"${manifest}\\"" >> "${status_file_name}"
    echo "    transport_parallel: \\"${transport_parallel}\\"" >> "${status_file_name}"
    echo "    transport_memory: \\"${transport_mem}\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"${exit_on_error_str}\\"" >> "${status_file_name}"

    # Add error message to status file if upload failed
    if [ \${UPLOAD_EXIT_CODE} -ne 0 ] && [ -n "\${ERROR_DETAILS:-}" ]; then
        echo "    error_details: |" >> "${status_file_name}"
        # Use printf to properly handle special characters and preserve formatting
        echo "\$ERROR_DETAILS" | while IFS= read -r line; do
            echo "        \$line" >> "${status_file_name}"
        done
    elif [ \${UPLOAD_EXIT_CODE} -ne 0 ]; then
        echo "    error_details: \\"SCORE upload failed\\"" >> "${status_file_name}"
    fi

    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        score-client: ${VERSION}
    END_VERSIONS

    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" = "true" ] && [ \${UPLOAD_EXIT_CODE} -ne 0 ]; then
        echo "SCORE upload failed and exit_on_error is true, exiting with error"
        exit \${UPLOAD_EXIT_CODE}
    else
        echo "Continuing workflow regardless of upload result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"
    """
    # Create output directory
    mkdir -p out
    
    # Create stub analysis_id file (copy input to output)
    cp ${analysis_id_file} ./out/analysis_id.txt
    
    # Create stub status file
    echo "process: \\"${task.process}\\"" > "${status_file_name}"
    echo "status: \\"SUCCESS\\"" >> "${status_file_name}"
    echo "exit_code: 0" >> "${status_file_name}"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${status_file_name}"
    echo "details:" >> "${status_file_name}"
    echo "    analysis_id: \\"stub-analysis-id-12345\\"" >> "${status_file_name}"
    echo "    file_manager_url: \\"${params.file_manager_url_upload ?: params.file_manager_url}\\"" >> "${status_file_name}"
    echo "    file_transfer_url: \\"${params.file_transfer_url_upload ?: params.file_transfer_url}\\"" >> "${status_file_name}"
    echo "    manifest_file: \\"${manifest}\\"" >> "${status_file_name}"
    echo "    transport_parallel: \\"${params.transport_parallel ?: task.cpus}\\"" >> "${status_file_name}"
    echo "    transport_memory: \\"${params.transport_mem ?: "2"}\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"false\\"" >> "${status_file_name}"

    # Create stub versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        score-client: stub
    END_VERSIONS
    """
}