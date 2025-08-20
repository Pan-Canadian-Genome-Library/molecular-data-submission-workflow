process CHECK_DEPENDENCIES {
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
        path "relational_mapping.json", emit: relational_mapping
        tuple (path("analysis_types.json")), emit : analysis_types
        path "versions.yml" , emit: versions

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
        ${workflow_file}

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
