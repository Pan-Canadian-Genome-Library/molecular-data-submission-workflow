process CHECK_CLINICAL {
    tag "${study_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:2.0.1--pyhdfd78af_0':
        'biocontainers/ampcombi:2.0.1--pyhdfd78af_0' }"

    input:
        tuple val(study_id)
        tuple val(token)
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
        --file_metadata ${file_metadata} \
        --analysis_metadata ${analysis_metadata} \
        --clinical_url ${params.clinical_url} \
        --file_manager_url ${params.file_manager_url} \
        --token ${token} \
        --study_id  ${study_id} \
        --data-directory ${data_directory} \
        --output-directory ${study_id} \
        ${sample_file} \
        ${specimen_file} \
        ${experiment_file} \
        ${read_group_file} \
        ${workflow_file}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkclinical: \$(python --version |& sed '1!d ; s/Python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${study_id}"
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkclinical: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
