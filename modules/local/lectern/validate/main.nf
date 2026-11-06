
process LECTERN_VALIDATE {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    //container 'ghcr.io/overture-stack/lectern'
    //container 'ghcr.io/pan-canadian-genome-library/clinical-submission:0.3.0'
    container 'ghcr.io/pan-canadian-genome-library/molecular-data-submission-workflow/lectern-validator:1.0.0'
    //containerOptions '--entrypoint ""'
    input:
    tuple val(meta), val(analysis), val(clinical), val(files), path(status_file), path(relational_mapping), path(analysis_types), path(data_directory)

    output:
    tuple val(meta), val(analysis), val(clinical), val(files) , path("*lectern_validate_status.yml") , emit : status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error.toString()
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def specimen_file = clinical.specimen!=null ? "${clinical.specimen}" : ""
    def sample_file = clinical.sample!=null ? "${clinical.sample}" : ""
    def experiment_file = clinical.experiment!=null ? "${clinical.experiment}" : ""
    def read_group_file = clinical.read_group!=null ? "${clinical.read_group}" : ""
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass' }" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping payload generation"
        GENERATION_EXIT_CODE=1
        ERROR_DETAILS="Submitting of clinical data skipped due to upstream failure"

    else
        echo "Upstream process successful, proceeding with payload generation"
        
        # Run main.py once and capture both exit code and error output
        app.ts \
        --url ${params.dictionary_manager_url} \
        --dictionary '${meta.study} pcgl schema' \
        --tsv \
        ${specimen_file} \
        ${sample_file} \
        ${experiment_file} \
        ${read_group_file} \
        --molecular \
        1> validation_results.log \
        2> generation_errors.tmp

        GENERATION_EXIT_CODE=\$?

        # Capture error details if generation failed
        ERROR_DETAILS=""
        if [ \$GENERATION_EXIT_CODE -ne 0 ]; then
            # Read all error details from captured stderr
            if [ -f "generation_errors.tmp" ] && [ -s "generation_errors.tmp" ]; then
                ERROR_DETAILS=\$(cat "generation_errors.tmp" | egrep 'ValueError\\: *.' || cat "generation_errors.tmp")
            else
                ERROR_DETAILS="Script execution failed"
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
        sample_meta : "${clinical.sample}"
        specimen_meta :  "${clinical.specimen}"
        experiment_meta : "${clinical.experiment}"
        read_group_meta : "${clinical.read_group}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message to status file if generation failed
    if [ \$GENERATION_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "    error_details: |" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
        echo "\$ERROR_DETAILS" | while IFS= read -r line; do
            echo "        \$line" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
        done
    fi
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        node: \$(node --version)
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
    
    // """
    // touch versions.yml
    // touch ${meta.id}.yml
    // #node bin/app.ts \
    // node --version > HELLO.txt
    // app.ts \
    // --url https://dictionary-manager.submission.genomelibrary.ca \
    // --dictionary 'PCGLST0002 pcgl schema' \
    // --tsv \
    // ${specimen_file} \
    // ${sample_file} \
    // ${experiment_file} \
    // ${read_group_file} \
    // 2> generation_errors.tmp

    // """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args
    
    touch ${prefix}.clinical_submission_status.yml 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lectern: \$(lectern --version)
    END_VERSIONS
    """
}
