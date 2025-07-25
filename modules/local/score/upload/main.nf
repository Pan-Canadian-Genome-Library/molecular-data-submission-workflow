
process SCORE_UPLOAD {
    tag "$meta.id"
    label 'process_medium'

    pod = [secret: workflow.runName + "-secret", mountPath: "/tmp/rdpc_secret"]
    container "${ params.score_container ?: 'ghcr.io/pan-canadian-genome-library/file-transfer' }:${ params.score_container_version ?: 'edge' }"

    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/score-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/score-client/logs"
    }

    input:
    tuple val(meta), path(analysis_id_file), path(manifest), path(upload)

    output:
    tuple val(meta), path("analysis_id.txt"),    emit: ready_to_publish
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    def song_url = params.song_url_upload ?: params.song_url
    def score_url = params.score_url_upload ?: params.score_url
    def transport_parallel = params.transport_parallel ?: task.cpus
    def transport_mem = params.transport_mem ?: "2"
    def accessToken = task.ext.api_upload_token ?: "`cat /tmp/rdpc_secret/secret`"
    def VERSION = params.score_container_version ?: '5.10.1'
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Read analysis ID from file
    ANALYSIS_ID=\$(cat ${analysis_id_file} | tr -d '\\n')
    
    # Copy analysis_id file to output for downstream processes
    cp ${analysis_id_file} analysis_id.txt
    
    # Check for upstream failures
    if [[ "\${ANALYSIS_ID}" == "ERROR:"* ]]; then
        echo "Detected failed upstream SONG submission, skipping SCORE upload"
        UPLOAD_EXIT_CODE=1
        ERROR_DETAILS="Skipped SCORE upload due to upstream SONG submission failure: \${ANALYSIS_ID}"
    elif grep -q "ERROR:" ${manifest} 2>/dev/null; then
        echo "Detected failed upstream manifest generation, skipping SCORE upload"
        UPLOAD_EXIT_CODE=1
        ERROR_DETAILS="Skipped SCORE upload due to upstream manifest generation failure"
    else
        echo "Valid inputs detected, proceeding with SCORE upload"
        
        export METADATA_URL=${song_url}
        export STORAGE_URL=${score_url}
        export TRANSPORT_PARALLEL=${transport_parallel}
        export TRANSPORT_MEMORY=${transport_mem}
        export ACCESSTOKEN=${accessToken}
        
        # Execute SCORE upload
        score-client upload --manifest ${manifest} $args
        UPLOAD_EXIT_CODE=\${?}
        
        if [ \${UPLOAD_EXIT_CODE} -ne 0 ]; then
            ERROR_DETAILS="SCORE upload command failed"
        fi
    fi

    # Create step-specific status file
    if [ \${UPLOAD_EXIT_CODE} -eq 0 ]; then
        STATUS_RESULT="SUCCESS"
    else
        STATUS_RESULT="FAILED"
    fi
    
    echo "process: \\"${task.process}\\"" > "${meta.id}_status.yml"
    echo "status: \\"\${STATUS_RESULT}\\"" >> "${meta.id}_status.yml"
    echo "exit_code: \${UPLOAD_EXIT_CODE}" >> "${meta.id}_status.yml"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${meta.id}_status.yml"
    echo "details:" >> "${meta.id}_status.yml"
    echo "    analysis_id: \\"\${ANALYSIS_ID}\\"" >> "${meta.id}_status.yml"
    echo "    song_url: \\"${song_url}\\"" >> "${meta.id}_status.yml"
    echo "    score_url: \\"${score_url}\\"" >> "${meta.id}_status.yml"
    echo "    manifest_file: \\"${manifest}\\"" >> "${meta.id}_status.yml"
    echo "    transport_parallel: \\"${transport_parallel}\\"" >> "${meta.id}_status.yml"
    echo "    transport_memory: \\"${transport_mem}\\"" >> "${meta.id}_status.yml"
    echo "    exit_on_error_enabled: \\"${exit_on_error_str}\\"" >> "${meta.id}_status.yml"
    
    # Add error message to status file if upload failed
    if [ \${UPLOAD_EXIT_CODE} -ne 0 ] && [ -n "\${ERROR_DETAILS:-}" ]; then
        echo "    error_message: \\"\${ERROR_DETAILS}\\"" >> "${meta.id}_status.yml"
    elif [ \${UPLOAD_EXIT_CODE} -ne 0 ]; then
        echo "    error_message: \\"SCORE upload failed\\"" >> "${meta.id}_status.yml"
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
}