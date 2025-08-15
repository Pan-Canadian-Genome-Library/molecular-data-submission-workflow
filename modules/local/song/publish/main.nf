process SONG_PUBLISH {
    tag "$meta.id"
    label 'process_single'

    container "${ params.file_manager_container ?: 'ghcr.io/pan-canadian-genome-library/file-manager-client' }:${ params.file_manager_container_tag ?: 'latest' }"
    containerOptions "-v \$(pwd):/song-client/logs"
  

    input:
    tuple val(meta), path(analysis_id_file)

    output:
    tuple val(meta), path("out/analysis_id.txt")             , emit: analysis_id
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    def file_manager_url = params.file_manager_url_upload ?: params.file_manager_url
    def accessToken = params.token
    def study_id = "${meta.study}"
    def VERSION = params.song_container_version ?: 'edge'
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"

    """
    # Set error handling to continue on failure for resilient processing
    set +e

    # Create output directory first
    mkdir -p out

    # Copy analysis_id file to output for downstream processes
    cp ${analysis_id_file} ./out/analysis_id.txt
    
    # Initialize variables
    PUBLISH_EXIT_CODE=0
    ERROR_DETAILS=""
    ANALYSIS_ID=""
    
    # Read analysis ID from file first (needed for both success and failure cases)
    if [ -f "${analysis_id_file}" ]; then
        ANALYSIS_ID=\$(cat ${analysis_id_file} | tr -d '\\n')
    fi
    
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping SONG publish"
        PUBLISH_EXIT_CODE=1
        ERROR_DETAILS="Skipped SONG publish due to upstream failure"
    else
        echo "Upstream process successful, proceeding with SONG publish"
        
        export CLIENT_SERVER_URL=${file_manager_url}
        export CLIENT_STUDY_ID=${study_id}
        export CLIENT_ACCESS_TOKEN=${accessToken}

        # Execute SONG publish
        sing publish -a \${ANALYSIS_ID} $args
        PUBLISH_EXIT_CODE=\${?}
        
        if [ \${PUBLISH_EXIT_CODE} -ne 0 ]; then
            # Capture error details preserving original formatting
            if [ -f ".command.err" ] && [ -s ".command.err" ]; then
                ERROR_DETAILS=\$(cat ".command.err")
            else
                ERROR_DETAILS="Script execution failed - no error details available"
            fi
        fi
    fi

    # Create step-specific status file
    if [ \${PUBLISH_EXIT_CODE} -eq 0 ]; then
        STATUS_RESULT="SUCCESS"
    else
        STATUS_RESULT="FAILED"
    fi
    
    echo "process: \\"${task.process}\\"" > "${status_file_name}"
    echo "status: \\"\${STATUS_RESULT}\\"" >> "${status_file_name}"
    echo "exit_code: \${PUBLISH_EXIT_CODE}" >> "${status_file_name}"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${status_file_name}"
    echo "details:" >> "${status_file_name}"
    echo "    study_id: \\"${study_id}\\"" >> "${status_file_name}"
    echo "    analysis_id: \\"\${ANALYSIS_ID}\\"" >> "${status_file_name}"
    echo "    file_manager_url: \\"${file_manager_url}\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"${exit_on_error_str}\\"" >> "${status_file_name}"

    # Add error message to status file if publish failed
    if [ \${PUBLISH_EXIT_CODE} -ne 0 ] && [ -n "\${ERROR_DETAILS:-}" ]; then
        echo "    error_details: |" >> "${status_file_name}"
        # Use printf to properly handle special characters and preserve formatting
        printf '%s\\n' "\$ERROR_DETAILS" | sed 's/^/        /' >> "${status_file_name}"
    elif [ \${PUBLISH_EXIT_CODE} -ne 0 ]; then
        echo "    error_details: \\"SONG publish failed\\"" >> "${status_file_name}"
    fi

    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: ${VERSION}
    END_VERSIONS

    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" = "true" ] && [ \${PUBLISH_EXIT_CODE} -ne 0 ]; then
        echo "SONG publish failed and exit_on_error is true, exiting with error"
        exit \${PUBLISH_EXIT_CODE}
    else
        echo "Continuing workflow regardless of publish result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"
    """
    # Create output directory first
    mkdir -p out
    # Create stub analysis_id file (copy input to output)
    cp ${analysis_id_file} ./out/analysis_id.txt

    # Create stub status file
    echo "process: \\"${task.process}\\"" > "${status_file_name}"
    echo "status: \\"SUCCESS\\"" >> "${status_file_name}"
    echo "exit_code: 0" >> "${status_file_name}"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${status_file_name}"
    echo "details:" >> "${status_file_name}"
    echo "    study_id: \\"${meta.study}\\"" >> "${status_file_name}"
    echo "    analysis_id: \\"stub-analysis-id-12345\\"" >> "${status_file_name}"
    echo "    file_manager_url: \\"${params.file_manager_url_upload ?: params.file_manager_url}\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"false\\"" >> "${status_file_name}"

    # Create stub versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: stub
    END_VERSIONS
    """
}