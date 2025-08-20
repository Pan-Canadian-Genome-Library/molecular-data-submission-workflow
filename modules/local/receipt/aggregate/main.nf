process RECEIPT_AGGREGATE {
    tag "$meta.batch_id"
    label 'process_single'

    conda "conda-forge::pyyaml=6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyyaml:6.0--py39h89e85a6_2' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(individual_receipts)

    output:
    tuple val(meta), path("*_batch_receipt.tsv"), emit: tsv_receipt
    tuple val(meta), path("*_batch_receipt.json"), emit: json_receipt
    tuple val(meta), path("*_batch_receipt.json"), path("*_batch_receipt.tsv"), emit: receipts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def timestamp = new Date().format("yyyyMMdd_HHmmss")
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

    # Generate batch receipt files using the external script
    main.py \\
        --individual-receipts ${individual_receipts.join(' ')} \\
        --batch-id "${timestamp}" \\
        --output-tsv "${timestamp}_batch_receipt.tsv" \\
        --output-json "${timestamp}_batch_receipt.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyyaml: \$(python -c "import yaml; print(yaml.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def timestamp = new Date().format("yyyyMMdd_HHmmss")
    """
    # Create mock batch receipt files
    cat <<-END_TSV > "${timestamp}_batch_receipt.tsv"
submitter_analysis_id	song_analysis_id	process	status	exit_code	timestamp	error_message	details	analysis_type	study_id	analysis_state	published_at
SAMPLE001	2fa0c2ab-b042-4ec4-a0c2-abb042bec48c	PAYLOAD_GENERATE	SUCCESS	0	2025-07-28T10:30:00Z		Payload generation completed	sequencing_experiment	TEST_STUDY	PUBLISHED	2025-07-28T16:12:59.392731
SAMPLE001	2fa0c2ab-b042-4ec4-a0c2-abb042bec48c	PAYLOAD_VALIDATE	SUCCESS	0	2025-07-28T10:31:00Z		Payload validation successful	sequencing_experiment	TEST_STUDY	PUBLISHED	2025-07-28T16:12:59.392731
SAMPLE002	3ba5e603-007f-6c0d-b471-48209f30d0fd	PAYLOAD_GENERATE	FAILED	1	2025-07-28T10:30:00Z	Missing required metadata fields	Payload generation failed	sequencing_experiment	TEST_STUDY	UNPUBLISHED	unknown
    END_TSV

    cat <<-END_JSON > "${timestamp}_batch_receipt.json"
{
  "batch_id": "${timestamp}",
  "generated_at": "2025-07-28T10:35:00Z",
  "total_analyses": 2,
  "successful_analyses": 1,
  "failed_analyses": 1,
  "analyses": [
    {
      "submitter_analysis_id": "SAMPLE001",
      "song_analysis_id": "2fa0c2ab-b042-4ec4-a0c2-abb042bec48c",
      "overall_status": "SUCCESS",
      "analysis_type": "sequencing_experiment",
      "study_id": "TEST_STUDY",
      "analysis_state": "PUBLISHED",
      "published_at": "2025-07-28T16:12:59.392731",
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
        }
      ]
    },
    {
      "submitter_analysis_id": "SAMPLE002",
      "song_analysis_id": "3ba5e603-007f-6c0d-b471-48209f30d0fd",
      "overall_status": "FAILED",
      "analysis_type": "sequencing_experiment",
      "study_id": "TEST_STUDY",
      "analysis_state": "UNPUBLISHED",
      "published_at": "unknown",
      "processes": [
        {
          "process": "PAYLOAD_GENERATE",
          "status": "FAILED",
          "exit_code": 1,
          "timestamp": "2025-07-28T10:30:00Z",
          "error_message": "Missing required metadata fields",
          "details": "Payload generation failed"
        }
      ]
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