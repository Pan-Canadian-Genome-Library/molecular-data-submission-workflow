//
// Subworkflow: File integrity validation with dedicated processes for each file type
// Processes files in parallel for maximum throughput
//

include { SAMTOOLS_QUICKCHECK      } from '../../../modules/local/samtools/quickcheck/main'
include { SEQKIT_SEQ             } from '../../../modules/local/seqkit/seq/main'
include { BCFTOOLS_VIEW          } from '../../../modules/local/bcftools/view/main'

workflow FILE_INTEGRITY {

    take:
    ch_payload      // channel: [val(meta), payload, [files]]

    main:
    ch_versions = Channel.empty()
    ch_status = Channel.empty()
    ch_validated_files = Channel.empty()

    // Flatten files for individual parallel processing
    // Each file will be processed separately for maximum parallelization  
    // Map from [meta, payload, [files]] to [meta.id, meta, payload, individual_file]
    ch_individual_files = ch_payload
        .flatMap { meta, payload, files ->
            // Handle case where files might be a single file or list of files
            def fileList = files instanceof List ? files : [files]
            fileList.collect { file -> [meta.id, meta, payload, file] }
        }
    
    // Separate individual files by type for dedicated validation processes
    // Group files with their associated index files for validation tools that need them
    
    // Store original file lists for index file lookup
    ch_original_files = ch_payload
        .map { meta, _payload, files ->
            [meta.id, files instanceof List ? files.flatten() : [files]]
        }
    
    // BAM/CRAM files with their index files (.bai, .csi)
    ch_bam_cram_files = ch_individual_files
        .filter { _meta_id, _meta, _payload, file -> 
            def name = file.name.toLowerCase()
            name.endsWith('.bam') || name.endsWith('.cram')
        }
        .combine(ch_original_files, by: 0)
        .map { _meta_id, meta, payload, main_file, all_files ->
            // Find associated index files for this BAM/CRAM file
            def index_files = all_files.findAll { index_file ->
                def index_name = index_file.name.toLowerCase()
                def main_name = main_file.name.toLowerCase()
                // Look for .bai, .csi, .crai index files that match the main file
                (index_name.startsWith(main_name) && (index_name.endsWith('.bai') || index_name.endsWith('.csi') || index_name.endsWith('.crai'))) ||
                (index_name == main_name + '.bai') || (index_name == main_name + '.csi') || (index_name == main_name + '.crai')
            }
            
            [meta, payload, main_file, index_files]
        }

    // VCF files with their index files (.tbi, .csi)
    ch_vcf_files = ch_individual_files
        .filter { _meta_id, _meta, _payload, file -> 
            def name = file.name.toLowerCase()
            (name.endsWith('.vcf') || name.endsWith('.vcf.gz') || name.endsWith('.bcf')) && 
            !(name.endsWith('.tbi') || name.endsWith('.csi'))
        }
        .combine(ch_original_files, by: 0)
        .map { _meta_id, meta, payload, main_file, all_files ->
            // Find associated index files for this VCF file
            def index_files = all_files.findAll { index_file ->
                def index_name = index_file.name.toLowerCase()
                def main_name = main_file.name.toLowerCase()
                // Look for .tbi, .csi index files that match the main file
                (index_name.startsWith(main_name) && (index_name.endsWith('.tbi') || index_name.endsWith('.csi'))) ||
                (index_name == main_name + '.tbi') || (index_name == main_name + '.csi')
            }
            
            [meta, payload, main_file, index_files]
        }

    // FASTQ files (no index files)
    ch_fastq_files = ch_individual_files
        .filter { _meta_id, _meta, _payload, file -> 
            def name = file.name.toLowerCase()
            name.endsWith('.fastq') || name.endsWith('.fastq.gz') || 
            name.endsWith('.fq') || name.endsWith('.fq.gz')
        }
        .map { _meta_id, meta, payload, file -> [meta, payload, file] }

    // Identify files that don't need validation (other files not covered above)
    // These will be passed through without validation
    // Note: Index files (.bai, .csi, .tbi) are now handled with their main files
    ch_passthrough_files = ch_individual_files
        .filter { _meta_id, _meta, _payload, file -> 
            def name = file.name.toLowerCase()
            // Files that are NOT BAM/CRAM/VCF/FASTQ and NOT index files (since index files are handled with main files)
            !(name.endsWith('.bam') || 
              name.endsWith('.cram') ||
              (name.endsWith('.vcf') || name.endsWith('.vcf.gz') || name.endsWith('.bcf')) ||
              name.endsWith('.fastq') || name.endsWith('.fastq.gz') || 
              name.endsWith('.fq') || name.endsWith('.fq.gz') ||
              name.endsWith('.bai') || name.endsWith('.csi') || name.endsWith('.tbi') || name.endsWith('.crai'))
        }
        .map { _meta_id, meta, payload, file -> [meta, payload, file] }

    // Debug file type separation
    ch_passthrough_files.view { meta, _payload, file ->
        "Passthrough file (no validation): meta=${meta.id}, file=${file.name}"
    }
    
    ch_bam_cram_files.view { meta, _payload, main_file, index_files ->
        "BAM/CRAM for validation: meta=${meta.id}, main=${main_file.name}, indexes=[${index_files.collect{it.name}.join(', ')}] (${index_files.size()} total)"
    }
    
    ch_vcf_files.view { meta, _payload, main_file, index_files ->
        "VCF for validation: meta=${meta.id}, main=${main_file.name}, indexes=[${index_files.collect{it.name}.join(', ')}] (${index_files.size()} total)"
    }

    // Run validation processes for individual files in parallel
    // Nextflow processes will handle empty channels gracefully
    // Use single samtools quickcheck process for both BAM and CRAM files
    SAMTOOLS_QUICKCHECK ( ch_bam_cram_files )
    ch_validated_files = ch_validated_files.mix(SAMTOOLS_QUICKCHECK.out.ch_validated_file)
    ch_status = ch_status.mix(SAMTOOLS_QUICKCHECK.out.status)
    ch_versions = ch_versions.mix(SAMTOOLS_QUICKCHECK.out.versions)

    BCFTOOLS_VIEW ( ch_vcf_files )
    ch_validated_files = ch_validated_files.mix(BCFTOOLS_VIEW.out.ch_validated_file)
    ch_status = ch_status.mix(BCFTOOLS_VIEW.out.status)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    SEQKIT_SEQ ( ch_fastq_files )
    ch_validated_files = ch_validated_files.mix(SEQKIT_SEQ.out.ch_validated_file)
    ch_status = ch_status.mix(SEQKIT_SEQ.out.status)
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    // Combine validated files with passthrough files (other files)
    // All files need to be regrouped together for downstream processing
    // Handle different output structures from validation modules:
    // - BAM/CRAM: [meta, payload, main_file, index_files[]]
    // - VCF: [meta, payload, main_file, index_files[]] 
    // - FASTQ: [meta, payload, file]
    // - Passthrough: [meta, payload, file]
    
    // Flatten validated files back to individual files
    ch_validated_flattened = ch_validated_files
        .flatMap { output ->
            if (output.size() == 4) {
                // BAM/CRAM/VCF with index files: [meta, payload, main_file, index_files[]]
                def meta = output[0]
                def payload = output[1]
                def main_file = output[2]
                def index_files = output[3] ?: []  // Handle null or empty index files
                
                // Start with main file, then add index files if they exist
                def all_files = [main_file]
                
                // Handle index_files whether it's a list, single file, or empty
                if (index_files != null && index_files != []) {
                    if (index_files instanceof List) {
                        // It's already a list - add all non-null files
                        index_files.findAll { it != null }.each { all_files.add(it) }
                    } else {
                        // It's a single file - add it
                        all_files.add(index_files)
                    }
                }

                // Return all files as individual tuples
                all_files.collect { file -> [meta, payload, file] }
            } else {
                // FASTQ: [meta, payload, file]
                [output]
            }
        }

    ch_all_files = ch_validated_flattened.mix(ch_passthrough_files)

    // Group all files (validated + passthrough) back by sample/meta for downstream processing
    // All files are now in format: [meta, payload, individual_file]
    // Need to collect all files back to: [meta, payload, [files]]
    ch_grouped_files = ch_all_files
        .map { meta, payload, file -> [meta.id, meta, payload, file] }
        .groupTuple(by: 0)
        .map { _id, metas, payloads, files ->
            // Use first meta and payload (should be identical for same sample)
            def meta = metas[0]
            def payload = payloads[0] 
            // Handle case where files might be nested arrays
            def flatFiles = files instanceof List ? files.flatten() : [files]
            [meta, payload, flatFiles]
        }

    // Debug final regrouped results
    ch_grouped_files.view { meta, payload, files ->
        "File_integrity: Final regrouped: meta=${meta.id}, payload=${payload.name}, files=[${files.collect{it.name}.join(', ')}] (${files.size()} total files)"
    }
    ch_status.view { meta, status ->
        "File_integrity: Status for ${meta.id}: ${status}"
    }

    emit:
    ch_payload_files = ch_grouped_files     // channel: [val(meta), payload, [files]]
    status           = ch_status            // channel: [val(meta), path(status)]
    versions         = ch_versions          // channel: [path(versions)]
}