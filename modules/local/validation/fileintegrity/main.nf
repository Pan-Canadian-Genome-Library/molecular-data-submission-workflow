process VALIDATION_FILEINTEGRITY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(payload), path(payload_files)

    output:
    tuple val(meta), path(payload), path(payload_files), emit: ch_payload_files
    tuple val(meta), path("*_status.yml"), emit: status
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"  // Convert boolean to string
    """
    # Set error handling to continue on failure for resilient processing
    set +e
    
    # Initialize variables
    INTEGRITY_EXIT_CODE=0
    ERROR_DETAILS=""
    
    # Check if upstream process was successful by checking meta.status
    if [ "${meta.status ?: 'pass'}" != "pass" ]; then
        echo "Upstream process failed (meta.status: ${meta.status ?: 'pass'}), skipping file integrity validation"
        INTEGRITY_EXIT_CODE=1
        ERROR_DETAILS="Skipped file integrity validation due to upstream failure"
    elif grep -q '"error".*"payload_generation_failed"' "${payload}" 2>/dev/null; then
        echo "Detected placeholder payload file, skipping file integrity validation"
        INTEGRITY_EXIT_CODE=1
        ERROR_DETAILS="Skipped file integrity validation - placeholder payload file detected"
    else
        echo "Upstream process successful, proceeding with file integrity validation"
        
        # Validate each file based on its type
        for file in ${payload_files}; do
            echo "Validating file: \$file"
            
            case "\$file" in
                *.bam)
                    echo "Validating BAM file: \$file"
                samtools quickcheck "\$file"
                if [ \$? -ne 0 ]; then
                    INTEGRITY_EXIT_CODE=1
                    ERROR_DETAILS="\${ERROR_DETAILS}BAM file validation failed for \$file; "
                fi
                ;;
            *.cram)
                echo "Validating CRAM file: \$file"
                samtools quickcheck "\$file"
                if [ \$? -ne 0 ]; then
                    INTEGRITY_EXIT_CODE=1
                    ERROR_DETAILS="\${ERROR_DETAILS}CRAM file validation failed for \$file; "
                fi
                ;;
            *.vcf|*.vcf.gz)
                echo "Validating VCF file: \$file"
                # Basic VCF header check
                if [[ "\$file" == *.gz ]]; then
                    header=\$(zcat "\$file" | head -1)
                else
                    header=\$(head -1 "\$file")
                fi
                if [[ "\$header" != "##fileformat=VCF"* ]]; then
                    INTEGRITY_EXIT_CODE=1
                    ERROR_DETAILS="\${ERROR_DETAILS}VCF file validation failed for \$file; "
                fi
                ;;
            *.fastq|*.fq|*.fastq.gz|*.fq.gz)
                echo "Validating FASTQ file: \$file"
                # Basic FASTQ format check - should start with @
                if [[ "\$file" == *.gz ]]; then
                    first_char=\$(zcat "\$file" | head -1 | cut -c1)
                else
                    first_char=\$(head -1 "\$file" | cut -c1)
                fi
                if [[ "\$first_char" != "@" ]]; then
                    INTEGRITY_EXIT_CODE=1
                    ERROR_DETAILS="\${ERROR_DETAILS}FASTQ file validation failed for \$file; "
                fi
                ;;
            *)
                echo "Warning: Unknown file type for \$file, skipping validation"
                ;;
        esac
    done
    fi
    
    # Create step-specific status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "\$(if [ \$INTEGRITY_EXIT_CODE -eq 0 ]; then echo 'SUCCESS'; else echo 'FAILED'; fi)"
    exit_code: \$INTEGRITY_EXIT_CODE
    timestamp: "\$(date -Iseconds)"
    details:
        analysis_id: "${meta.id}"
        files_validated: "${payload_files}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS
    
    # Add error message to status file if validation failed
    if [ \$INTEGRITY_EXIT_CODE -ne 0 ] && [ -n "\$ERROR_DETAILS" ]; then
        echo "        error_message: \"\$ERROR_DETAILS\"" >> "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    fi
    
    # Always create versions.yml before any exit
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    
    # Only exit with error if exit_on_error is explicitly true
    if [ "${exit_on_error_str}" = "true" ] && [ \$INTEGRITY_EXIT_CODE -ne 0 ]; then
        echo "File integrity validation failed and exit_on_error is true, exiting with error"
        exit \$INTEGRITY_EXIT_CODE
    else
        echo "Continuing workflow regardless of integrity validation result (exit_on_error=${exit_on_error_str})"
        exit 0
    fi
    """

    stub:
    def exit_on_error = params.exit_on_error ?: task.ext.exit_on_error ?: false
    def exit_on_error_str = exit_on_error ? "true" : "false"
    """
    # Create mock integrity status file
    cat <<-END_STATUS > "${meta.id}_${task.process.toLowerCase().replace(':', '_')}_status.yml"
    process: "${task.process}"
    status: "SUCCESS"
    exit_code: 0
    timestamp: "2025-01-25T10:30:00+00:00"
    details:
        analysis_id: "${meta.id}"
        files_validated: "${payload_files}"
        exit_on_error_enabled: "${exit_on_error_str}"
    END_STATUS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.17
    END_VERSIONS
    """
}
