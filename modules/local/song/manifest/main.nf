process SONG_MANIFEST {
    tag "$meta.id"
    label 'process_single'

    container "${ params.song_container ?: 'ghcr.io/pan-canadian-genome-library/file-manager-client' }:${ params.song_container_version ?: 'edge' }"

    if (workflow.containerEngine == "singularity") {
        containerOptions "--bind \$(pwd):/song-client/logs"
    } else if (workflow.containerEngine == "docker") {
        containerOptions "-v \$(pwd):/song-client/logs"
    }

    input:
    tuple val(meta), path(analysis_id_file), path(payload), path(upload)

    output:
    path "out/manifest.txt"       , emit: manifest
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml"           , emit: versions
    tuple val(meta), path("out/analysis_id.txt"), path("out/manifest.txt"), path(upload),    emit: manifest_upload

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    def song_url = params.song_url_upload ?: params.song_url
    def accessToken = params.token
    def VERSION = params.song_container_version ?: 'edge'
    def study_id = "${meta.study_id}"
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Create output directory first
    mkdir -p out
    
    # Initialize variables
    MANIFEST_EXIT_CODE=0
    ERROR_DETAILS=""
    ANALYSIS_ID=""
    
    # Read analysis ID from file first (needed for both success and failure cases)
    if [ -f "${analysis_id_file}" ]; then
        ANALYSIS_ID=\$(cat ${analysis_id_file} | tr -d '\\n')
    fi
    
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping manifest generation"
        MANIFEST_EXIT_CODE=1
        ERROR_DETAILS="Skipped manifest generation due to upstream failure"
        
        # Create placeholder manifest file
        echo "# Manifest generation skipped due to upstream failure" > out/manifest.txt
    else
        echo "Upstream process successful, proceeding with manifest generation"
        
        export CLIENT_SERVER_URL=${song_url}
        export CLIENT_STUDY_ID=${study_id}
        export CLIENT_ACCESS_TOKEN=${accessToken}

        # Execute SONG manifest generation
        sing manifest -a \${ANALYSIS_ID} -d . -f out/manifest.txt $args
        MANIFEST_EXIT_CODE=\${?}
        
        if [ \${MANIFEST_EXIT_CODE} -ne 0 ]; then
            ERROR_DETAILS="SONG manifest command failed"
            # Create placeholder manifest file for failed generation
            echo "# Manifest generation failed" > out/manifest.txt
        fi
    fi

    # Copy analysis_id file to output for downstream processes
    cp ${analysis_id_file} ./out/analysis_id.txt

    # Create step-specific status file
    if [ \${MANIFEST_EXIT_CODE} -eq 0 ]; then
        STATUS_RESULT="SUCCESS"
    else
        STATUS_RESULT="FAILED"
    fi
    
    echo "process: \\"${task.process}\\"" > "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "status: \\"\${STATUS_RESULT}\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "exit_code: \${MANIFEST_EXIT_CODE}" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "details:" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "    study_id: \\"${study_id}\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "    analysis_id: \\"\${ANALYSIS_ID}\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "    song_url: \\"${song_url}\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "    manifest_file: \\"./out/manifest.txt\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    echo "    exit_on_error_enabled: \\"${exit_on_error_str}\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    
    # Add error message to status file if manifest generation failed
    if [ \${MANIFEST_EXIT_CODE} -ne 0 ] && [ -n "\${ERROR_DETAILS:-}" ]; then
        echo "    error_message: \\"\${ERROR_DETAILS}\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    elif [ \${MANIFEST_EXIT_CODE} -ne 0 ]; then
        echo "    error_message: \\"Manifest generation failed\\"" >> "${meta.id}_${task.process.toLowerCase()}_status.yml"
    fi

    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: ${VERSION}
    END_VERSIONS

    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" = "true" ] && [ \${MANIFEST_EXIT_CODE} -ne 0 ]; then
        echo "Manifest generation failed and exit_on_error is true, exiting with error"
        exit \${MANIFEST_EXIT_CODE}
    else
        echo "Continuing workflow regardless of manifest result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """
}