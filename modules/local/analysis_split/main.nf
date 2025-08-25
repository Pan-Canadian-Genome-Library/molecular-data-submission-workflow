process ANALYSIS_SPLIT {
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
        tuple path(relational_mapping) //Json

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    path("*/*/*status.yml"), emit: status
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${study_id}"
    def sample_file = sample_metadata ? "--sample_metadata ${sample_metadata}" : ""
    def specimen_file = specimen_metadata ? "--specimen_metadata ${specimen_metadata}" : ""
    def experiment_file = experiment_metadata ? "--experiment_metadata ${experiment_metadata}" : ""
    def read_group_file = read_group_metadata ? "--read_group_metadata ${read_group_metadata}" : ""
    def workflow_file = workflow_metadata ? "--workflow_metadata ${workflow_metadata}" : ""
    """
    main.py \
        --relational_mapping ${relational_mapping} \
        --file_metadata ${file_metadata} \
        --analysis_metadata ${analysis_metadata} \
        ${workflow_file} \
        ${sample_file} \
        ${specimen_file} \
        ${experiment_file} \
        ${read_group_file} \
        --output_directory ${study_id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${study_id}"
    """
    mkdir -p ${study_id}/{EXAMPLE_01,EXAMPLE_02,EXAMPLE_03}

    touch ${study_id}/EXAMPLE_01/{analysis.tsv,files.tsv,workflow.tsv,specimen.tsv,sample.tsv,experiment.tsv,read_group.tsv}
    touch ${study_id}/EXAMPLE_02/{analysis.tsv,files.tsv,specimen.tsv,sample.tsv,experiment.tsv}
    touch ${study_id}/EXAMPLE_03/{analysis.tsv,files.tsv}

    cat <<-END_STATUS > "${study_id}/EXAMPLE_01/EXAMPLE_01_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-22T10:30:00+00:00"
    details:
    details:
        analysis_id: "${study_id}}" 
        sample_meta : "${sample_metadata}"
        specimen_meta :  "${specimen_metadata}"
        experiment_meta : "${experiment_metadata}"
        read_group_meta : "${read_group_metadata}"
    END_STATUS

    cat <<-END_STATUS > "${study_id}/EXAMPLE_02/EXAMPLE_02_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-22T10:30:00+00:00"
    details:
    details:
        analysis_id: "${study_id}}" 
        sample_meta : "${sample_metadata}"
        specimen_meta :  "${specimen_metadata}"
        experiment_meta : "${experiment_metadata}"
        read_group_meta : "${read_group_metadata}"
    END_STATUS

    cat <<-END_STATUS > "${study_id}/EXAMPLE_03/EXAMPLE_03_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-22T10:30:00+00:00"
    details:
    details:
        analysis_id: "${study_id}}" 
        sample_meta : "${sample_metadata}"
        specimen_meta :  "${specimen_metadata}"
        experiment_meta : "${experiment_metadata}"
        read_group_meta : "${read_group_metadata}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        analysissplit: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
