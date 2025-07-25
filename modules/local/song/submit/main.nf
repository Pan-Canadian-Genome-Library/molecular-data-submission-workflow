
process SONG_SUBMIT {
    tag "$meta.id"
    label 'process_single'

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]

    container "${ params.song_container ?: 'ghcr.io/pan-canadian-genome-library/file-manager-client' }:${ params.song_container_version ?: 'edge' }"
    
    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/song-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/song-client/logs"
    }

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
    def accessToken = task.ext.api_upload_token ?: "`cat /tmp/rdpc_secret/secret`"
    def VERSION = params.song_container_version ?: 'edge'
    def study_id = "${meta.study_id}"
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    export CLIENT_SERVER_URL=${song_url}
    export CLIENT_STUDY_ID=${study_id}
    export CLIENT_ACCESS_TOKEN=${accessToken}

    # Execute SONG submit and capture output
    ANALYSIS_ID=\$(sing submit -f ${payload} $args | jq -er .analysisId | tr -d '\\n' 2>&1)
    SUBMIT_EXIT_CODE=\${?}
    
    # Create placeholder analysis ID if submission failed
    if [ \${SUBMIT_EXIT_CODE} -ne 0 ]; then
        ANALYSIS_ID="ERROR: SONG submission failed"
    fi
    
    # Write ANALYSIS_ID to file for Nextflow to capture
    echo "\${ANALYSIS_ID}" > analysis_id.txt

    # Create step-specific status file
    if [ \${SUBMIT_EXIT_CODE} -eq 0 ]; then
        STATUS_RESULT="SUCCESS"
    else
        STATUS_RESULT="FAILED"
    fi
    
    echo "process: \\"${task.process}\\"" > "${meta.id}_status.yml"
    echo "status: \\"\${STATUS_RESULT}\\"" >> "${meta.id}_status.yml"
    echo "exit_code: \${SUBMIT_EXIT_CODE}" >> "${meta.id}_status.yml"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${meta.id}_status.yml"
    echo "details:" >> "${meta.id}_status.yml"
    echo "    study_id: \\"${study_id}\\"" >> "${meta.id}_status.yml"
    echo "    analysis_id: \\"\${ANALYSIS_ID}\\"" >> "${meta.id}_status.yml"
    echo "    song_url: \\"${song_url}\\"" >> "${meta.id}_status.yml"
    echo "    payload_file: \\"${payload}\\"" >> "${meta.id}_status.yml"
    echo "    exit_on_error_enabled: \\"${exit_on_error_str}\\"" >> "${meta.id}_status.yml"
    
    # Add error message to status file if submission failed
    if [ \${SUBMIT_EXIT_CODE} -ne 0 ]; then
        echo "    error_message: \\"\${ANALYSIS_ID}\\"" >> "${meta.id}_status.yml"
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