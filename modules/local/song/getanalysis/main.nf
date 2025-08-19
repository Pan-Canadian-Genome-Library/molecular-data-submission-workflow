// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process SONG_GETANALYSIS {
    tag "$meta.id"
    label 'process_single'

    container "${ params.file_manager_container ?: 'ghcr.io/pan-canadian-genome-library/file-manager-client' }:${ params.file_manager_container_tag ?: 'latest' }"

    // Container options must be set as a string, not a Groovy closure
    containerOptions "-v \$(pwd):/song-client/logs"

    input:
    tuple val(meta), path(analysis_id_file)

    output:
    tuple val(meta), path("*.analysis.json"), emit: analysis_json
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    def file_manager_url = params.file_manager_url_download ?: params.file_manager_url
    def accessToken = params.token
    def VERSION = params.song_container_version ?: 'edge'
    def study_id = "${meta.study}"
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Initialize variables
    GETANALYSIS_EXIT_CODE=0
    ERROR_DETAILS=""
    ANALYSIS_ID=""
    
    # Read analysis ID from file first (needed for both success and failure cases)
    if [ -f "${analysis_id_file}" ]; then
        ANALYSIS_ID=\$(cat ${analysis_id_file} | tr -d '\\n')
    fi
    
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping analysis retrieval"
        GETANALYSIS_EXIT_CODE=1
        ERROR_DETAILS="Skipped analysis retrieval due to upstream failure"
        
        # Create placeholder analysis file with ANALYSIS_ID prefix
        echo '{"error": "Analysis retrieval skipped due to upstream failure"}' > ${meta.id}_\${ANALYSIS_ID}.analysis.json
    else
        echo "Upstream process successful, proceeding with analysis retrieval"
        
        export CLIENT_SERVER_URL=${file_manager_url}
        export CLIENT_STUDY_ID=${study_id}
        export CLIENT_ACCESS_TOKEN=${accessToken}

        # Execute SONG get analysis with ANALYSIS_ID prefix
        sing search -a \${ANALYSIS_ID} $args > ${meta.id}_\${ANALYSIS_ID}.analysis.json
        GETANALYSIS_EXIT_CODE=\${?}
        
        if [ \${GETANALYSIS_EXIT_CODE} -ne 0 ]; then
            # Capture error details preserving original formatting
            if [ -f ".command.err" ] && [ -s ".command.err" ]; then
                ERROR_DETAILS=\$(cat ".command.err")
            else
                ERROR_DETAILS="Script execution failed - no error details available"
            fi

            # Create placeholder analysis file for failed retrieval with ANALYSIS_ID prefix
            echo '{"error": "Analysis retrieval failed"}' > ${meta.id}_\${ANALYSIS_ID}.analysis.json
        fi
    fi

    # Create step-specific status file
    if [ \${GETANALYSIS_EXIT_CODE} -eq 0 ]; then
        STATUS_RESULT="SUCCESS"
    else
        STATUS_RESULT="FAILED"
    fi
    
    echo "process: \\"${task.process}\\"" > "${status_file_name}"
    echo "status: \\"\${STATUS_RESULT}\\"" >> "${status_file_name}"
    echo "exit_code: \${GETANALYSIS_EXIT_CODE}" >> "${status_file_name}"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${status_file_name}"
    echo "details:" >> "${status_file_name}"
    echo "    study_id: \\"${study_id}\\"" >> "${status_file_name}"
    echo "    analysis_id: \\"\${ANALYSIS_ID}\\"" >> "${status_file_name}"
    echo "    file_manager_url: \\"${file_manager_url}\\"" >> "${status_file_name}"
    echo "    analysis_file: \\"\${ANALYSIS_ID}.analysis.json\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"${exit_on_error_str}\\"" >> "${status_file_name}"

    # Add error message to status file if retrieval failed
    if [ \${GETANALYSIS_EXIT_CODE} -ne 0 ] && [ -n "\${ERROR_DETAILS:-}" ]; then
        echo "    error_details: |" >> "${status_file_name}"
        echo "\$ERROR_DETAILS" | sed 's/^/            /' >> "${status_file_name}"
    elif [ \${GETANALYSIS_EXIT_CODE} -ne 0 ]; then
        echo "    error_details: \\"Analysis retrieval failed\\"" >> "${status_file_name}"
    fi

    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: ${VERSION}
    END_VERSIONS

    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" = "true" ] && [ \${GETANALYSIS_EXIT_CODE} -ne 0 ]; then
        echo "Analysis retrieval failed and exit_on_error is true, exiting with error"
        exit \${GETANALYSIS_EXIT_CODE}
    else
        echo "Continuing workflow regardless of retrieval result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def status_file_name = "${meta.id}_" + (task.process.toLowerCase().replace(':', '_')) + "_status.yml"
    """
    # Create stub analysis file
    echo '{"analysisId": "stub-analysis-id-12345", "studyId": "stub-study", "state": "PUBLISHED"}' > stub-analysis-id-12345.analysis.json
    
    # Create stub status file
    echo "process: \\"${task.process}\\"" > "${status_file_name}"
    echo "status: \\"SUCCESS\\"" >> "${status_file_name}"
    echo "exit_code: 0" >> "${status_file_name}"
    echo "timestamp: \\"\$(date -Iseconds)\\"" >> "${status_file_name}"
    echo "details:" >> "${status_file_name}"
    echo "    study_id: \\"${meta.study}\\"" >> "${status_file_name}"
    echo "    analysis_id: \\"${meta.analysis_id}\\"" >> "${status_file_name}"
    echo "    file_manager_url: \\"${params.file_manager_url_upload ?: params.file_manager_url}\\"" >> "${status_file_name}"
    echo "    exit_on_error_enabled: \\"false\\"" >> "${status_file_name}"

    # Create stub versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        song-client: stub
    END_VERSIONS
    """
}
