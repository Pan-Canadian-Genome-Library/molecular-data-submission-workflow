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
    tuple val({ meta.findAll { key, value -> key != 'status' } }), path("*_status.yml"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = task.ext.exit_on_error ?: false
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    generate_payload.py \\
        --submitter-analysis-id "${meta.id}" \\
        --analysis-type "${meta.type}" \\
        --study-id "${meta.study}" \\
        --file-meta "${file_meta}" \\
        --analysis-meta "${analysis_meta}" \\
        --workflow-meta "${workflow_meta}" \\
        --data-files ${data_files} \\
        --output "${prefix}_payload.json"

    GENERATION_EXIT_CODE=\$?

    # Capture stderr for error details if generation failed
    ERROR_DETAILS=""
    if [ \$GENERATION_EXIT_CODE -ne 0 ]; then
        ERROR_DETAILS=\$(generate_payload.py \\
            --submitter-analysis-id "${meta.id}" \\
            --analysis-type "${meta.type}" \\
            --study-id "${meta.study}" \\
            --file-meta "${file_meta}" \\
            --analysis-meta "${analysis_meta}" \\
            --workflow-meta "${workflow_meta}" \\
            --data-files ${data_files} \\
            --output "${prefix}_payload.json" 2>&1 | tail -1 | xargs || echo "Script execution failed")
            
        # Create empty/placeholder payload file if generation failed to satisfy Nextflow output requirements
        if [ ! -f "${prefix}_payload.json" ]; then
            echo '{"error": "payload_generation_failed", "message": "Placeholder file created due to generation failure"}' > "${prefix}_payload.json"
        fi
    fi

    # Create step-specific status file (meta.id used for grouping in channel)
    cat <<-END_STATUS > "${meta.id}_status.yml"
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
        exit_on_error_enabled: "${exit_on_error}"
    END_STATUS
    
    # Add error message to status file if generation failed
    if [ \$GENERATION_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "        error_message: \"\$ERROR_DETAILS\"" >> "${meta.id}_status.yml"
    fi
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    
    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error}" == "true" ] && [ \$GENERATION_EXIT_CODE -ne 0 ]; then
        echo "Payload generation failed and exit_on_error is true, exiting with error"
        exit \$GENERATION_EXIT_CODE
    else
        echo "Continuing workflow regardless of generation result (exit_on_error=${exit_on_error})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = task.ext.exit_on_error ?: false
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create mock TSV metadata files for testing with correct field structures
    
    # File metadata
    echo -e "submitter_analysis_id\\tfileName\\tfileSize\\tfileMd5sum\\tfileType\\tfileAccess\\tdataType" > mock_file_meta.tsv
    echo -e "${meta.id}\\ttest_file.fastq.gz\\t1000000\\tabc123def456\\tFASTQ\\tcontrolled\\tAligned Reads" >> mock_file_meta.tsv
    
    # Analysis metadata  
    echo -e "studyId\\tsubmitter_analysis_id\\tanalysisType\\tsubmitter_participant_id\\tsubmitter_specimen_id\\tsubmitter_sample_id\\tsubmitter_experiment_id\\tdata_category\\tvariant_class\\tvariant_calling_strategy\\tgenome_build\\tgenome_annotation" > mock_analysis_meta.tsv
    echo -e "${meta.study}\\t${meta.id}\\t${meta.type}\\tPART_001\\tSPEC_001\\tSAMPLE_001\\tEXP_001\\tgenomic\\tSNV\\tgatk\\tGRCh38\\tGENCODE_v39" >> mock_analysis_meta.tsv
    
    # Workflow metadata
    echo -e "submitter_workflow_id\\tsubmitter_analysis_id\\tworkflow_name\\tworkflow_version\\tworkflow_url" > mock_workflow_meta.tsv
    echo -e "WF_001\\t${meta.id}\\ttest_workflow\\t1.0.0\\thttps://example.com/workflow" >> mock_workflow_meta.tsv
    
    generate_payload.py \\
        --submitter-analysis-id "${meta.id}" \\
        --analysis-type "${meta.type}" \\
        --study-id "${meta.study}" \\
        --file-meta mock_file_meta.tsv \\
        --analysis-meta mock_analysis_meta.tsv \\
        --workflow-meta mock_workflow_meta.tsv \\
        --output "${prefix}_payload.json"

    # Create mock step-specific status file
    cat <<-END_STATUS > "${meta.id}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-22T10:30:00+00:00"
    details:
        analysis_id: "${meta.id}"
        payload_file: "${prefix}_payload.json"
        file_meta: "mock_file_meta.tsv"
        analysis_meta: "mock_analysis_meta.tsv"
        workflow_meta: "mock_workflow_meta.tsv"
        exit_on_error_enabled: "${exit_on_error}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
