

process CLINICALSUBMISSION {
    tag "$meta.id"
    label 'process_single'

container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:2.0.1--pyhdfd78af_0':
        'biocontainers/ampcombi:2.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(analysis), val(clinical), val(files), path(status_file)
    val token
    val clinical_url

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.bam"), emit: bam
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error.toString()
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_file = clinical.sample ? "--sample_metadata ${clinical.sample}" : ""
    def specimen_file = clinical.experiment ? "--specimen_metadata ${clinical.experiment}" : ""
    def experiment_file = clinical.experiment ? "--experiment_metadata ${clinical.experiment}" : ""
    def read_group_file = clinical.read_group ? "--read_group_metadata ${clinical.read_group}" : ""
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
        main.py \\
        --sample_metadata ${sample_file} \
        --specimen_metadata  ${specimen_file} \
        --experiment_metadata ${experiment_file} \
        --read_group_metadata ${read_group_file} \
        --clinical_url ${clinical_url} \
        --token ${token} \
        --study_id ${meta.study_id} 2>generation_errors.tmp

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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clinicalsubmission: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
