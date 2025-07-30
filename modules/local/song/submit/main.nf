
process SONG_SUBMIT {
    tag "$meta.id"
    label 'process_single'

    container "${ params.song_container ?: 'ghcr.io/pan-canadian-genome-library/file-manager-client' }:${ params.song_container_version ?: 'edge' }"
    containerOptions "-v \$(pwd):/song-client/logs"


    input:
    tuple val(meta), path(payload), path(files)

    output:
    tuple val(meta), path("analysis_id.txt"), emit: analysis_id
    tuple val(meta), path("analysis_id.txt"), path(payload), path(files), emit: analysis_files
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    def song_url = params.song_url_upload ?: params.song_url
    def accessToken = params.token
    def VERSION = params.song_container_version ?: 'edge'
    def study_id = "${meta.study}"
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"

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
        
        export CLIENT_SERVER_URL=${song_url}
        export CLIENT_STUDY_ID=${study_id}
        export CLIENT_ACCESS_TOKEN=${accessToken}

        # Execute SONG submit and capture output
        ANALYSIS_ID=\$(sing submit -f ${payload} $args | jq -er .analysisId | tr -d '\\n' 2>&1)
        SUBMIT_EXIT_CODE=\${?}
        
        # Create placeholder analysis ID if submission failed
        if [ \${SUBMIT_EXIT_CODE} -ne 0 ]; then
            ANALYSIS_ID="unavailable"
            ERROR_DETAILS="SONG submit command failed"
        fi
    fi
    
    # Write ANALYSIS_ID to file for Nextflow to capture
    echo "\${ANALYSIS_ID}" > analysis_id.txt

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
    echo "    song_url: \\"${song_url}\\"" >> "${status_file_name}"
    echo "    payload_file: \\"${payload}\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"${exit_on_error_str}\\"" >> "${status_file_name}"

    # Add error message to status file if submission failed
    if [ \${SUBMIT_EXIT_CODE} -ne 0 ] && [ -n "\${ERROR_DETAILS:-}" ]; then
        echo "    error_message: \\"\${ERROR_DETAILS}\\"" >> "${status_file_name}"
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
}