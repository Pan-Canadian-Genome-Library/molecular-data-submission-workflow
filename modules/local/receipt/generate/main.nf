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

process RECEIPT_GENERATE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pyyaml=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyyaml:6.0--py39h89e85a6_2' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(status_files), path(analysis_file, optional: true)

    output:
    tuple val(meta), path("${meta.id}_receipt.json"), emit: json_receipt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install required Python packages to temporary directory
    echo "Installing required Python packages..."
    
    # Create a temporary directory for package installation
    TEMP_PYTHON_LIB="\$(mktemp -d)/python_packages"
    mkdir -p "\$TEMP_PYTHON_LIB"
    
    # Install to temporary directory and add to Python path
    pip install --target "\$TEMP_PYTHON_LIB" PyYAML
    if [ \$? -ne 0 ]; then
        echo "Failed to install PyYAML package" >&2
        exit 1
    fi
    
    # Set PYTHONPATH to include our temporary package directory
    export PYTHONPATH="\$TEMP_PYTHON_LIB:\${PYTHONPATH:-}"


    # Generate receipt files using the external script
    META_ARGS="--submitter-analysis-id \"${meta.id}\" --study-id \"${meta.study}\" --analysis-type \"${meta.type}\""
    if [[ -f "${analysis_file}" ]]; then
        main.py \
            --status-files ${status_files.join(' ')} \
            --analysis-file "${analysis_file}" \
            $META_ARGS \
            --output-json "${meta.id}_receipt.json"
    else
        main.py \
            --status-files ${status_files.join(' ')} \
            $META_ARGS \
            --output-json "${meta.id}_receipt.json"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    # Create mock receipt file
    cat <<-END_JSON > "${meta.id}_receipt.json"
{
  "submitter_analysis_id": "SAMPLE001",
  "song_analysis_id": "2fa0c2ab-b042-4ec4-a0c2-abb042bec48c",
  "overall_status": "SUCCESS",
  "analysis_type": "sequencing_experiment",
  "study_id": "TEST_STUDY",
  "analysis_state": "PUBLISHED",
  "published_at": "2025-07-28T16:12:59.392731",
  "generated_at": "2025-07-28T10:35:00Z",
  "processes": [
    {
      "process": "PAYLOAD_GENERATE",
      "status": "SUCCESS",
      "exit_code": 0,
      "timestamp": "2025-07-28T10:30:00Z",
      "details": "Payload generation completed"
    },
    {
      "process": "PAYLOAD_VALIDATE", 
      "status": "SUCCESS",
      "exit_code": 0,
      "timestamp": "2025-07-28T10:31:00Z",
      "details": "Payload validation successful"
    },
    {
      "process": "SONG_SUBMIT",
      "status": "SUCCESS", 
      "exit_code": 0,
      "timestamp": "2025-07-28T10:32:00Z",
      "details": "Analysis submitted to SONG"
    }
  ]
}
    END_JSON

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        pyyaml: "6.0"
    END_VERSIONS
    """
}
