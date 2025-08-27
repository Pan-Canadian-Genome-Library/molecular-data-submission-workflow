
process SONG_SUBMIT {
    tag "$meta.id"
    label 'process_single'

    container "${params.file_manager_container}:${params.file_manager_container_tag}"
    containerOptions "-v \$(pwd):/song-client/logs"


    input:
    tuple val(meta), path(payload), path(files)

    output:
    tuple val(meta), path("*_analysis_id.txt"), emit: analysis_id
    tuple val(meta), path("*_analysis_id.txt"), path(payload), path(files), emit: analysis_files
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    def file_manager_url = params.file_manager_url
    def accessToken = params.token
    def VERSION = params.file_manager_container_tag
    def study_id = "${meta.study}"
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"
    def allow_duplicates_arg  = params.allow_duplicates ? "--allow-duplicates" : ''
    def analysis_id_file = "${meta.id}_analysis_id.txt"
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Initialize variables
    SUBMIT_EXIT_CODE=0
    ERROR_DETAILS=""
    
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping SONG submit"
        SUBMIT_EXIT_CODE=1
        ERROR_DETAILS="Skipped SONG submit due to upstream failure"
        ANALYSIS_ID="unavailable"
    else
        echo "Upstream process successful, proceeding with SONG submit"

        export CLIENT_SERVER_URL=${file_manager_url}
        export CLIENT_STUDY_ID=${study_id}
        export CLIENT_ACCESS_TOKEN=${accessToken}

        # Execute SONG submit and capture output
        ANALYSIS_ID=\$(sing submit -f ${payload} $args $allow_duplicates_arg | jq -er .analysisId | tr -d '\\n' 2>&1)
        SUBMIT_EXIT_CODE=\${?}
        
        # Create placeholder analysis ID if submission failed
        if [ \${SUBMIT_EXIT_CODE} -ne 0 ]; then
            ANALYSIS_ID="unavailable"
            # Capture error details preserving original formatting
            if [ -f "song.log" ] && [ -s "song.log" ]; then
                ERROR_DETAILS=\$(cat "song.log")
            else
                ERROR_DETAILS="Script execution failed - no error details available"
            fi
        fi
    fi
    
    # Write ANALYSIS_ID to file for Nextflow to capture
    echo "\${ANALYSIS_ID}" > "${analysis_id_file}"

    # Create step-specific status file
    if [ \${SUBMIT_EXIT_CODE} -eq 0 ]; then
        STATUS_RESULT="SUCCESS"
    else
        STATUS_RESULT="FAILED"
    fi

    echo "process: \\"${task.process}\\"" > "${status_file_name}"
    echo "status: \\"\${STATUS_RESULT}\\"" >> "${status_file_name}"
    echo "exit_code: \${SUBMIT_EXIT_CODE}" >> "${status_file_name}"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${status_file_name}"
    echo "details:" >> "${status_file_name}"
    echo "    study_id: \\"${study_id}\\"" >> "${status_file_name}"
    echo "    analysis_id: \\"\${ANALYSIS_ID}\\"" >> "${status_file_name}"
    echo "    file_manager_url: \\"${file_manager_url}\\"" >> "${status_file_name}"
    echo "    payload_file: \\"${payload}\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"${exit_on_error_str}\\"" >> "${status_file_name}"

    # Add error message to status file if submission failed
    if [ \${SUBMIT_EXIT_CODE} -ne 0 ] && [ -n "\${ERROR_DETAILS:-}" ]; then
        echo "    error_details: |" >> "${status_file_name}"
        echo "\$ERROR_DETAILS" | sed 's/^/            /' >> "${status_file_name}"

    elif [ \${SUBMIT_EXIT_CODE} -ne 0 ]; then
        echo "    error_message: \\"Submission failed\\"" >> "${status_file_name}"
    fi

    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: ${VERSION}
    END_VERSIONS

    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" = "true" ] && [ \${SUBMIT_EXIT_CODE} -ne 0 ]; then
        echo "SONG submission failed and exit_on_error is true, exiting with error"
        exit \${SUBMIT_EXIT_CODE}
    else
        echo "Continuing workflow regardless of submission result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi

    """

    stub:
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"
    def analysis_id_file = "${meta.id}_analysis_id.txt"
    """
    # Create stub analysis ID
    echo "stub-analysis-id-12345" > "${analysis_id_file}"

    # Create stub status file
    echo "process: \\"${task.process}\\"" > "${status_file_name}"
    echo "status: \\"SUCCESS\\"" >> "${status_file_name}"
    echo "exit_code: 0" >> "${status_file_name}"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${status_file_name}"
    echo "details:" >> "${status_file_name}"
    echo "    study_id: \\"${meta.study}\\"" >> "${status_file_name}"
    echo "    analysis_id: \\"stub-analysis-id-12345\\"" >> "${status_file_name}"
    echo "    file_manager_url: \\"${params.file_manager_url_upload ?: params.file_manager_url}\\"" >> "${status_file_name}"
    echo "    payload_file: \\"${payload}\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"false\\"" >> "${status_file_name}"

    # Create stub versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: stub
    END_VERSIONS
    """
}