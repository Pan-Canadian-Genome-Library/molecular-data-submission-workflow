process CHECK_CLINICAL {
    tag "${study_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:2.0.1--pyhdfd78af_0':
        'biocontainers/ampcombi:2.0.1--pyhdfd78af_0' }"

    input:
        tuple val(study_id)
        tuple path(file_metadata) // Spreadsheet
        tuple path(analysis_metadata) // Spreadsheet
        tuple path(workflow_metadata) // Spreadsheet
        tuple path(read_group_metadata) // Spreadsheet
        tuple path(experiment_metadata) // Spreadsheet
        tuple path(specimen_metadata) // Spreadsheet
        tuple path(sample_metadata) // Spreadsheet
        tuple path(data_directory) // file path

    output:
        path "*/*/*.{SUCCESS,FAILURE}" , emit: status
        path "versions.yml" , emit: versions
        path "*_status.yml", emit : status_yml

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error.toString()
    def prefix = task.ext.prefix ?: "${study_id}"
    def args = task.ext.args ?: ''
    def sample_file = sample_metadata ? "--sample_metadata ${sample_metadata}" : ""
    def specimen_file = specimen_metadata ? "--specimen_metadata ${specimen_metadata}" : ""
    def experiment_file = experiment_metadata ? "--experiment_metadata ${experiment_metadata}" : ""
    def read_group_file = read_group_metadata ? "--read_group_metadata ${read_group_metadata}" : ""
    def workflow_file = workflow_metadata ? "--workflow_metadata ${workflow_metadata}" : ""
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Check if upstream process was successful by checking meta.status
    echo "Proceeding with dependency check"
    
    # Run main.py once and capture both exit code and error output
    main.py \
        --file_metadata ${file_metadata} \
        --analysis_metadata ${analysis_metadata} \
        --clinical_url ${params.clinical_url} \
        --file_manager_url ${params.file_manager_url} \
        --token ${params.token} \
        --study_id  ${study_id} \
        --data-directory ${params.data_directory} \
        --output-directory ${study_id} \
        ${sample_file} \
        ${specimen_file} \
        ${experiment_file} \
        ${read_group_file} \
        ${workflow_file} 2>generation_errors.tmp

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

    # Create step-specific status file (meta.id used for grouping in channel)
    cat <<-END_STATUS > "${study_id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$GENERATION_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$GENERATION_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    details:
        analysis_id: "${study_id}"
        payload_file: "${prefix}_payload.json"
        file_meta: "${file_metadata}"
        analysis_meta: "${analysis_metadata}"
        workflow_meta: "${workflow_metadata}"
        sample_meta: "${sample_metadata}"
        specimen_meta: "${specimen_metadata}"
        experiment_meta: "${experiment_metadata}"
        read_group_meta: "${read_group_metadata}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message to status file if generation failed
    if [ \$GENERATION_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "    error_message: \\"\$ERROR_DETAILS\\"" >> "${study_id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
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
    def prefix = task.ext.prefix ?: "${study_id}"
    def args = task.ext.args ?: ''
    def sample_file = sample_metadata ? "--sample_metadata ${sample_metadata}" : ""
    def specimen_file = specimen_metadata ? "--specimen_metadata ${specimen_metadata}" : ""
    def experiment_file = experiment_metadata ? "--experiment_metadata ${experiment_metadata}" : ""
    def read_group_file = read_group_metadata ? "--read_group_metadata ${read_group_metadata}" : ""
    def workflow_file = workflow_metadata ? "--workflow_metadata ${workflow_metadata}" : ""
    """
        mkdir -p ${study_id}/example_analysis_01
        touch ${study_id}/example_analysis_01/registerable_analysis.tsv
        touch ${study_id}/example_analysis_01/registerable_files.tsv
        touch ${study_id}/example_analysis_01/registerable_workflow.tsv
        touch ${study_id}/example_analysis_01/registerable_sample.tsv
        touch ${study_id}/example_analysis_01/registerable_specimen.tsv
        touch ${study_id}/example_analysis_01/registerable_experiment.tsv
        touch ${study_id}/example_analysis_01/registerable_read_group.tsv
        touch ${study_id}/example_analysis_01/example_analysis_01.SUCCESS
        mkdir -p ${study_id}/example_analysis_02
        touch ${study_id}/example_analysis_02/unregisterable_analysis.tsv
        touch ${study_id}/example_analysis_02/unregisterable_files.tsv
        touch ${study_id}/example_analysis_02/unregisterable_workflow.tsv
        touch ${study_id}/example_analysis_02/unregisterable_sample.tsv
        touch ${study_id}/example_analysis_02/unregisterable_specimen.tsv
        touch ${study_id}/example_analysis_02/unregisterable_experiment.tsv
        touch ${study_id}/example_analysis_02/unregisterable_read_group.tsv
        touch ${study_id}/example_analysis_02/example_analysis_01.FAILURE

    # Create mock step-specific status file
    cat <<-END_STATUS > "${study_id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-22T10:30:00+00:00"
    details:
        analysis_id: "${study_id}"
        payload_file: "${prefix}_payload.json"
        file_meta: "${file_metadata}"
        analysis_meta: "${analysis_metadata}"
        workflow_meta: "${workflow_metadata}"
        sample_meta: "${sample_metadata}"
        specimen_meta: "${specimen_metadata}"
        experiment_meta: "${experiment_metadata}"
        read_group_meta: "${read_group_metadata}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
