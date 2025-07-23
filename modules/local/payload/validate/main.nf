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

process PAYLOAD_VALIDATE {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'ghcr.io/icgc-argo-workflows/data-processing-utility-tools.payload-gen-seq-experiment:0.8.3' }"

    input:
    tuple val(meta), path(payload_file), path(data_files)

    output:
    tuple val(meta), path(payload_file), path(data_files), emit: validated_payload_files
    tuple val({ meta.findAll { key, value -> key != 'status' } }), path("*_status.yml"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = task.ext.exit_on_error ?: false
    def schema_url = task.ext.schema_url ?: ''
    
    if (!schema_url) {
        error "schema_url must be provided via task.ext.schema_url"
    }
    
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Check if upstream payload generation was successful by checking meta.status
    if [ "${meta.status}" != "success" ]; then
        echo "Upstream process failed (meta.status: ${meta.status}), skipping validation"
        VALIDATION_EXIT_CODE=1
        ERROR_DETAILS="Skipped validation due to upstream failure"
    elif grep -q '"error": "payload_generation_failed"' "${payload_file}" 2>/dev/null; then
        echo "Detected placeholder payload file, skipping validation"
        VALIDATION_EXIT_CODE=1
        ERROR_DETAILS="Skipped validation - placeholder payload file detected"
    else
        echo "Upstream process successful, proceeding with validation"
        
        validate_payload.py \\
            --payload "${payload_file}" \\
            --schema-url "${schema_url}" \\
            --output-report "validation_report_${meta.id}.txt"
        
        VALIDATION_EXIT_CODE=\$?
        
        # Read validation report for error details if validation failed
        ERROR_DETAILS=""
        if [ \$VALIDATION_EXIT_CODE -ne 0 ] && [ -f "validation_report_${meta.id}.txt" ]; then
            # Extract only MESSAGE lines from the validation report
            ERROR_DETAILS=\$(grep "^MESSAGE:" "validation_report_${meta.id}.txt" | sed 's/^MESSAGE: //' | tr '\\n' ' | ' | sed 's/ | \$//')
            
            # If no MESSAGE lines found, fall back to last line
            if [ -z "\$ERROR_DETAILS" ]; then
                ERROR_DETAILS=\$(tail -1 "validation_report_${meta.id}.txt" | xargs)
            fi
        fi
    fi
    
    # Create step-specific status file (meta.id used for grouping in channel)
    cat <<-END_STATUS > "${meta.id}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$VALIDATION_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$VALIDATION_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    details:
        analysis_id: "${meta.id}"
        payload_file: "${payload_file}"
        schema_url: "${schema_url}"
        exit_on_error_enabled: "${exit_on_error}"
    END_STATUS
    
    # Add error message to status file if validation failed
    if [ \$VALIDATION_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "        error_message: \"\$ERROR_DETAILS\"" >> "${meta.id}_status.yml"
    fi
    
    # Update meta with analysis status for downstream decision making
    # This will be used in the Nextflow process to update the meta map
    echo "\$VALIDATION_EXIT_CODE" > validation_exit_code.txt
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    
    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error}" == "true" ] && [ \$VALIDATION_EXIT_CODE -ne 0 ]; then
        echo "Validation failed and exit_on_error is true, exiting with error"
        exit \$VALIDATION_EXIT_CODE
    else
        echo "Continuing workflow regardless of validation result (exit_on_error=${exit_on_error})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = task.ext.exit_on_error ?: false
    def schema_url = task.ext.schema_url ?: 'https://example.com/schema.json'
    """
    # Create mock step-specific status file
    cat <<-END_STATUS > "${meta.id}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-22T10:30:00+00:00"
    details:
        analysis_id: "${meta.id}"
        payload_file: "${payload_file}"
        schema_url: "${schema_url}"
        exit_on_error_enabled: "${exit_on_error}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
