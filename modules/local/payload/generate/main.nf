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

process PAYLOAD_GENERATE {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:    
    tuple val(meta), path(file_meta), path(analysis_meta), path(workflow_meta), path(data_files)
    // where file_meta, analysis_meta, workflow_meta are TSV files

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*_payload.json"), path(data_files), emit: payload_files  
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error.toString()
    def prefix = task.ext.prefix ?: "${meta.id}"
    def workflow_meta_arg = (workflow_meta && workflow_meta != [] && workflow_meta.toString() != '[]') ? "--workflow-meta \"${workflow_meta}\"" : ""
    """
    # Set error handling to continue on failure for resilient processing
    set +e

    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping payload generation"
        GENERATION_EXIT_CODE=1
        ERROR_DETAILS="Skipped payload generation due to upstream failure"

        # Create placeholder payload file to satisfy Nextflow output requirements
        echo '{"error": "payload_generation_failed", "message": "Placeholder file created due to upstream failure"}' > "${prefix}_payload.json"
    else
        echo "Upstream process successful, proceeding with payload generation"

        # Run main.py once and capture both exit code and error output
        main.py \
            --submitter-analysis-id "${meta.id}" \
            --analysis-type "${meta.type}" \
            --study-id "${meta.study}" \
            --file-meta "${file_meta}" \
            --analysis-meta "${analysis_meta}" \
            ${workflow_meta_arg} \
            --data-files ${data_files} \
            --output "${prefix}_payload.json" 2>generation_errors.tmp

        GENERATION_EXIT_CODE=\$?

        # Capture error details if generation failed
        ERROR_DETAILS=""
        if [ \$GENERATION_EXIT_CODE -ne 0 ]; then
            # Read all error details from captured stderr
            if [ -f "generation_errors.tmp" ] && [ -s "generation_errors.tmp" ]; then
                ERROR_DETAILS=\$(cat "generation_errors.tmp" | tr '\\n' ' | ' | sed 's/ | \$//' || echo "Script execution failed")
            else
                ERROR_DETAILS="Script execution failed"
            fi

            # Create empty/placeholder payload file if generation failed to satisfy Nextflow output requirements
            if [ ! -f "${prefix}_payload.json" ]; then
                echo '{"error": "payload_generation_failed", "message": "Placeholder file created due to generation failure"}' > "${prefix}_payload.json"
            fi
        fi
    fi

    # Create step-specific status file (meta.id used for grouping in channel)
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$GENERATION_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$GENERATION_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}"
        payload_file: "${prefix}_payload.json"
        file_meta: "${file_meta}"
        analysis_meta: "${analysis_meta}"
        workflow_meta: "${workflow_meta}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message to status file if generation failed
    if [ \$GENERATION_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "    error_message: \\"\$ERROR_DETAILS\\"" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    fi
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    
    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" == "true" ] && [ \$GENERATION_EXIT_CODE -ne 0 ]; then
        echo "Payload generation failed and exit_on_error is true, exiting with error"
        exit \$GENERATION_EXIT_CODE
    else
        echo "Continuing workflow regardless of generation result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error.toString()
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create simple mock payload JSON file
    echo '{"mock": "payload", "analysis_id": "${meta.id}", "study_id": "${meta.study}", "analysis_type": "${meta.type}"}' > "${prefix}_payload.json"

    # Create mock step-specific status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-22T10:30:00+00:00"
    work_directory: "\$PWD"
    details:
        analysis_id: "${meta.id}"
        payload_file: "${prefix}_payload.json"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
