#!/usr/bin/env python3
"""
Aggregate individual analysis receipts into batch receipts (JSON + TSV).
"""

import argparse
import json
import csv
from datetime import datetime
import sys


def load_individual_receipts(receipt_files):
    """Load and parse all individual receipt JSON files."""
    receipts = []
    
    for file_path in receipt_files:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                receipts.append(data)
        except Exception as e:
            print(f"Warning: Could not load receipt file {file_path}: {e}", file=sys.stderr)
    
    return receipts


def generate_batch_tsv(receipts, batch_id, output_file):
    """Generate batch TSV receipt from individual receipts."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Header
        writer.writerow([
            'submitter_analysis_id', 'overall_status', 'process', 'status', 'exit_code', 
            'timestamp', 'error_message', 
            'file_manager_analysis_id', 'analysis_type', 'study_id', 'analysis_state', 'published_at'
        ])
        
        # Data rows
        for receipt in receipts:
            submitter_analysis_id = receipt.get('submitter_analysis_id', None)
            file_manager_analysis_id = receipt.get('file_manager_analysis_id', None)
            overall_status = receipt.get('overall_status', None)
            analysis_type = receipt.get('analysis_type', None)
            study_id = receipt.get('study_id', None)
            analysis_state = receipt.get('analysis_state', None)
            published_at = receipt.get('published_at', None)
            process_failure_point = None

            if receipt.get("processes"):
                failed_processes=[]
                for process in receipt.get('processes'):
                    if process.get('status'):
                        if process.get('status')==False:
                            failed_processes.append(process.get("process"))
                
                if not overall_status:
                    process_failure_point = failed_processes[-1]

            for process in receipt.get('processes', []):
                writer.writerow([
                    submitter_analysis_id,
                    overall_status,
                    process.get('process', None),
                    process.get('status', None),
                    process.get('exit_code', None),
                    process.get('timestamp', None),
                    process.get('error_message', None),
                    file_manager_analysis_id,
                    analysis_type,
                    study_id,
                    analysis_state,
                    published_at
                ])


def generate_batch_json(receipts, batch_id, output_file):
    """Generate batch JSON receipt from individual receipts."""
    batch_receipt = {
        'batch_id': batch_id,
        'generated_at': datetime.now().isoformat() + 'Z',
        'total_analyses': len(receipts),
        'successful_analyses': 0,
        'failed_analyses': 0,
        'analyses': []
    }
    
    for receipt in receipts:
        # Count successful/failed analyses
        if receipt.get('overall_status') == 'SUCCESS':
            batch_receipt['successful_analyses'] += 1
        else:
            batch_receipt['failed_analyses'] += 1
        
        # Add the individual receipt to the batch
        batch_receipt['analyses'].append(receipt)
    
    with open(output_file, 'w') as f:
        json.dump(batch_receipt, f, indent=2)
    
    # Generate summary files for easy parsing by Nextflow
    with open('total.txt', 'w') as f:
        f.write(str(batch_receipt['total_analyses']))
    
    with open('success.txt', 'w') as f:
        f.write(str(batch_receipt['successful_analyses']))
        
    with open('failed.txt', 'w') as f:
        f.write(str(batch_receipt['failed_analyses']))


def main():
    parser = argparse.ArgumentParser(description='Aggregate individual receipts into batch receipt')
    parser.add_argument('--individual-receipts', nargs='+', required=True,
                        help='List of individual receipt JSON files')
    parser.add_argument('--batch-id', required=True,
                        help='Batch identifier (timestamp)')
    parser.add_argument('--output-tsv', required=True,
                        help='Output TSV file path')
    parser.add_argument('--output-json', required=True,
                        help='Output JSON file path')
    
    args = parser.parse_args()
    
    # Load all individual receipts
    print(f"Loading {len(args.individual_receipts)} individual receipts", file=sys.stderr)
    receipts = load_individual_receipts(args.individual_receipts)
    
    if not receipts:
        print("No valid receipt files found", file=sys.stderr)
        sys.exit(1)
    
    # Generate batch outputs
    generate_batch_tsv(receipts, args.batch_id, args.output_tsv)
    generate_batch_json(receipts, args.batch_id, args.output_json)
    
    print(f"Generated batch receipt for {len(receipts)} analyses", file=sys.stderr)
    print(f"Batch ID: {args.batch_id}", file=sys.stderr)
    print(f"TSV: {args.output_tsv}", file=sys.stderr)
    print(f"JSON: {args.output_json}", file=sys.stderr)


if __name__ == '__main__':
    main()
