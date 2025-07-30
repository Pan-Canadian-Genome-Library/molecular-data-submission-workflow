//
// Subworkflow for generating batch receipts from status files and analysis data
//
// Usage examples:
// 1. All analyses in one batch (current implementation):
//    BATCH_RECEIPT_GENERATION(status_analysis_ch)
//
// 2. To group by specific field, modify the grouping logic in the main section
//    For example, to group by study_id or batch_id from meta
//

include { RECEIPT_GENERATE } from '../../../modules/local/receipt/generate/main'
include { RECEIPT_AGGREGATE } from '../../../modules/local/receipt/aggregate/main'

workflow BATCH_RECEIPT_GENERATION {
    
    take:
    ch_status_analysis    // channel: tuple val(meta), path(status_files), path(analysis_file)
    
    main:
    
    // Generate individual receipts for each analysis
    RECEIPT_GENERATE (
        ch_status_analysis
    )
    
    // Option 1: Collect all receipts into a single batch (current implementation)
    all_receipts = RECEIPT_GENERATE.out.json_receipt
        .map { _meta, receipt_file -> receipt_file }
        .collect()
        .map { receipt_files ->
            // Create batch meta with timestamp-based batch_id
            def timestamp = new Date().format('yyyyMMdd_HHmmss')
            def batch_meta = [
                id: "batch_${timestamp}",
                batch_id: "batch_${timestamp}"
            ]
            return [batch_meta, receipt_files]
        }
    
    
    // Aggregate individual receipts into batch receipts
    RECEIPT_AGGREGATE (
        all_receipts
    )
    
    emit:
    batch_tsv_receipt   = RECEIPT_AGGREGATE.out.tsv_receipt    // tuple val(meta), path(batch_receipt.tsv)
    batch_json_receipt  = RECEIPT_AGGREGATE.out.json_receipt   // tuple val(meta), path(batch_receipt.json)
    versions            = RECEIPT_GENERATE.out.versions.mix(RECEIPT_AGGREGATE.out.versions)
    
}
