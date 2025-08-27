process ENHANCE_RECEIPT_METADATA {
    tag "$meta.batch_id"
    label 'process_single'

    input:
    tuple val(meta), path(batch_receipt_tsv), path(batch_receipt_json)
    path(batch_summary_json)

    output:
    tuple val(meta), path(batch_receipt_tsv), emit: tsv_receipt
    tuple val(meta), path(batch_receipt_json), emit: json_receipt
    tuple val(meta), path(batch_receipt_json), path(batch_receipt_tsv), emit: receipts

    script:
    """
    # Read the batch summary JSON and extract metadata
    python3 << 'EOF'
import json
import sys

# Read the batch summary file
try:
    with open('${batch_summary_json}', 'r') as f:
        summary = json.load(f)
    
    # Create enhanced metadata and write to a file
    enhanced_meta = {
        'batch_id': summary.get('batch_id', '${meta.batch_id ?: "unknown"}'),
        'study': summary.get('study', '${meta.study ?: "unknown"}'),
        'total_analyses': summary.get('total_analyses', 0),
        'successful_analyses': summary.get('successful_analyses', 0),
        'failed_analyses': summary.get('failed_analyses', 0),
        'success_rate': summary.get('success_rate', 0.0),
        'status': summary.get('status', 'UNKNOWN')
    }
    
    # Save enhanced metadata for use by Nextflow
    with open('enhanced_meta.json', 'w') as f:
        json.dump(enhanced_meta, f, indent=2)
        
    print(f"Enhanced metadata: {enhanced_meta}")
    
except Exception as e:
    print(f"Error reading batch summary: {e}", file=sys.stderr)
    # Create default metadata
    default_meta = {
        'batch_id': '${meta.batch_id ?: "unknown"}',
        'study': '${meta.study ?: "unknown"}',
        'total_analyses': 0,
        'successful_analyses': 0,
        'failed_analyses': 0,
        'success_rate': 0.0,
        'status': 'ERROR'
    }
    with open('enhanced_meta.json', 'w') as f:
        json.dump(default_meta, f, indent=2)
EOF
    """
}
