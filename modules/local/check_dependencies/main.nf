process CHECK_DEPENDENCIES {
    tag "${study_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:2.0.1--pyhdfd78af_0':
        'biocontainers/ampcombi:2.0.1--pyhdfd78af_0' }"

    input:
        val(study_id)
        path(file_metadata) // Spreadsheet
        path(analysis_metadata) // Spreadsheet
        path(workflow_metadata) // Spreadsheet
        path(read_group_metadata) // Spreadsheet
        path(experiment_metadata) // Spreadsheet
        path(specimen_metadata) // Spreadsheet
        path(sample_metadata) // Spreadsheet
        path(data_directory) // file path
    output:
        path "relational_mapping.json", emit: relational_mapping
        path "analysis_types.json", emit : analysis_types
        path "versions.yml" , emit: versions
        tuple val(study_id), path("*_status.yml"), emit: status

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${study_id}"
    def sample_file = sample_metadata ? "--sample_metadata ${sample_metadata}" : ""
    def specimen_file = specimen_metadata ? "--specimen_metadata ${specimen_metadata}" : ""
    def experiment_file = experiment_metadata ? "--experiment_metadata ${experiment_metadata}" : ""
    def read_group_file = read_group_metadata ? "--read_group_metadata ${read_group_metadata}" : ""
    def workflow_file = workflow_metadata ? "--workflow_metadata ${workflow_metadata}" : ""
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Initialize variables
    CHECKDEPENDENCIES_EXIT_CODE=0
    ERROR_DETAILS=""

    main.py \
        --file_metadata ${file_metadata} \
        --analysis_metadata ${analysis_metadata} \
        --clinical_url ${params.clinical_url} \
        --file_manager_url ${params.file_manager_url} \
        --study_id  ${study_id} \
        --token ${params.token} \
        ${sample_file} \
        ${specimen_file} \
        ${experiment_file} \
        ${read_group_file} \
        ${workflow_file} 2> check_dependencies_errors.tmp

    CHECKDEPENDENCIES_EXIT_CODE=\$?

    # Capture error details if the command failed
    if [ \$CHECKDEPENDENCIES_EXIT_CODE -ne 0 ]; then
        ERROR_DETAILS=\$(cat check_dependencies_errors.tmp)

        # Create empty output files if they don't exist due to failure
        [ ! -f "relational_mapping.json" ] && echo '{}' > relational_mapping.json
        [ ! -f "analysis_types.json" ] && echo '{}' > analysis_types.json
    fi
    
    # Create step-specific status file
    cat <<-END_STATUS > "${study_id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$CHECKDEPENDENCIES_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$CHECKDEPENDENCIES_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    work_directory: "\$PWD"
    details:
        study_id: "${study_id}"
        file_metadata: "${file_metadata}"
        analysis_metadata: "${analysis_metadata}"
        workflow_metadata: "${workflow_metadata}"
        read_group_metadata: "${read_group_metadata}"
        experiment_metadata: "${experiment_metadata}"
        specimen_metadata: "${specimen_metadata}"
        sample_metadata: "${sample_metadata}"
        files_directory: "${data_directory}"
    END_STATUS
    
    # Add error message to status file if check dependencies failed
    if [ \$CHECKDEPENDENCIES_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        # Format multi-line error details properly for YAML
        echo "    error_details: |" >> "${study_id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
        echo "\$ERROR_DETAILS" | sed 's/^/        /' >> "${study_id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    fi
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${study_id}"
    def sample_file = sample_metadata ? "--sample_metadata ${sample_metadata}" : ""
    def specimen_file = specimen_metadata ? "--specimen_metadata ${specimen_metadata}" : ""
    def experiment_file = experiment_metadata ? "--experiment_metadata ${experiment_metadata}" : ""
    def read_group_file = read_group_metadata ? "--read_group_metadata ${read_group_metadata}" : ""
    def workflow_file = workflow_metadata ? "--workflow_metadata ${workflow_metadata}" : ""
    """
    touch relational_mapping.json
    touch analysis_types.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
