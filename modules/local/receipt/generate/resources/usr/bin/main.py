#!/usr/bin/env python3
"""
Aggregate status information from multiple YAML status files and analysis into batch receipts.
Expects unique status filenames from different processes.
"""

import argparse
import yaml
import json
import csv
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import sys


def parse_status_file(file_path):
    """Parse a single status YAML file."""
    try:
        with open(file_path, 'r') as f:
            data = yaml.safe_load(f)
        
        # Extract analysis_id from various possible sources (this is submitter_analysis_id)
        submitter_analysis_id = None
        if 'details' in data and 'analysis_id' in data['details']:
            submitter_analysis_id = data['details']['analysis_id']
        elif 'analysis_id' in data:
            submitter_analysis_id = data['analysis_id']
        else:
            # Try to extract from filename if not in YAML
            filename = Path(file_path).stem
            if '_status' in filename:
                # Handle patterns like: SAMPLE001_metadata_status, SAMPLE001_crosscheck_status
                parts = filename.split('_')
                if len(parts) >= 2 and parts[-1] == 'status':
                    # Remove process name and '_status' suffix
                    submitter_analysis_id = '_'.join(parts[:-2]) if len(parts) > 2 else parts[0]
                else:
                    submitter_analysis_id = filename.replace('_status', '')
        
        return {
            'submitter_analysis_id': submitter_analysis_id,
            'process': data.get('process', None),
            'status': data.get('status', None),
            'exit_code': data.get('exit_code', None),
            'timestamp': data.get('timestamp', None),
            'work_directory': data.get('work_directory', None),
            'details': data.get('details', None),
            'source_file': str(file_path)
        }
    except Exception as e:
        print(f"Error parsing {file_path}: {e}", file=sys.stderr)
        return None


def simplfy_error_msg(og_msg):
    error_code_msg={
        "403" : "403 - Invalid token found!",
        "500" : "500 - Error on server end. Please contact PCGL Admin for assistance.",
        "bio.overture.song.sdk.errors.ServerResponseErrorHandler.handleError" : "Error on sever end (FileManager/FileTrasfer). Please contact PCGL Admin for assistance.",
        "401" : "401 - User is unauthorized to perform this action. Please contact PCGL Admin for assistance."
    }
    for line in og_msg.split("\n"):
        for key in error_code_msg.keys():
            if key in line:
                return(error_code_msg[key])

    return(og_msg)

def generate_individual_receipt(processes, analysis_info, output_file):
    """Generate individual analysis receipt in JSON format."""
    # Determine overall status for this analysis
    overall_status = 'SUCCESS'
    if any(p['status'] == 'FAILED' for p in processes):
        overall_status = 'FAILED'
    
    # Clean up process data for JSON output
    clean_processes = []
    for p in processes:

        if p.get('details'):
            if p.get('details')!=None:
                if p.get('details').get('error_details'):
                    if p.get('details').get('error_details')!=None:
                        p['details']['error_details']=simplfy_error_msg(p.get('details').get('error_details'))

        clean_process = {
            'process': p['process'],
            'status': p['status'],
            'exit_code': p['exit_code'],
            'timestamp': p['timestamp'],
            'work_directory': p.get('work_directory', None),
            'details':   p.get('details', None),
        }

        clean_processes.append(clean_process)
    
    # Create individual receipt structure
    receipt = {
        'submitter_analysis_id': analysis_info.get('submitter_analysis_id', processes[0]['submitter_analysis_id']),
        'file_manager_analysis_id': analysis_info.get('file_manager_analysis_id', None),
        'overall_status': overall_status,
        'process_failure_point':  [process.get('process') for process in clean_processes if process.get('status')=='FAILED'][0] if overall_status=='FAILED' else None,
        'analysis_type': analysis_info.get('analysis_type', None),
        'study_id': analysis_info.get('study_id', None),
        'analysis_state': analysis_info.get('file_manager_analysis_state', None),
        'published_at': analysis_info.get('file_manager_analysis_publish_at', None),
        'generated_at': datetime.now().isoformat() + 'Z',
        'processes': clean_processes
    }
    
    with open(output_file, 'w') as f:
        json.dump(receipt, f, indent=2)


def main():
    parser = argparse.ArgumentParser(description='Generate individual analysis receipt')
    parser.add_argument('--status-files', nargs='+', required=True,
                        help='List of status YAML files')
    parser.add_argument('--analysis-file', required=False,
                        help='Analysis JSON file (optional, only for file_manager_analysis_id, song_analysis_state, song_analysis_publish_at)')
    parser.add_argument('--submitter-analysis-id', required=True, help='Submitter analysis id (from meta)')
    parser.add_argument('--study-id', required=True, help='Study id (from meta)')
    parser.add_argument('--analysis-type', required=True, help='Analysis type (from meta)')
    parser.add_argument('--output-json', required=True,
                        help='Output JSON file path')

    args = parser.parse_args()

    # Parse analysis file if provided, else use defaults for song fields
    file_manager_analysis_id = 'Not applicable'
    file_manager_analysis_state = 'Not applicable'
    file_manager_analysis_publish_at = 'Not applicable'
    if args.analysis_file:
        try:
            with open(args.analysis_file, 'r') as f:
                data = json.load(f)
            file_manager_analysis_id = data.get('analysisId', None)
            file_manager_analysis_state = data.get('analysisState', None)
            file_manager_analysis_publish_at = data.get('publishedAt', None)
        except Exception as e:
            print(f"Error parsing analysis file {args.analysis_file}: {e}", file=sys.stderr)
    else:
        print("No analysis file provided, using default values for file_manager_analysis_id, file_manager_analysis_state, file_manager_analysis_publish_at.", file=sys.stderr)

    # Build analysis_info dict for receipt
    analysis_info = {
        'file_manager_analysis_id': file_manager_analysis_id,
        'file_manager_analysis_state': file_manager_analysis_state,
        'file_manager_analysis_publish_at': file_manager_analysis_publish_at,
        'submitter_analysis_id': args.submitter_analysis_id,
        'study_id': args.study_id,
        'analysis_type': args.analysis_type
    }
    # Parse all status files (they are already grouped by submitter_analysis_id)
    processes = []
    submitter_analysis_id = None

    print(f"Processing {len(args.status_files)} status files", file=sys.stderr)

    for status_file in args.status_files:
        data = parse_status_file(status_file)
        if data and data.get('submitter_analysis_id'):
            if submitter_analysis_id is None:
                submitter_analysis_id = data['submitter_analysis_id']
            elif submitter_analysis_id != data['submitter_analysis_id']:
                print(f"Warning: Mismatched submitter_analysis_ids detected: {submitter_analysis_id} vs {data['submitter_analysis_id']}", file=sys.stderr)
            processes.append(data)
        else:
            print(f"Warning: Could not extract submitter_analysis_id from {status_file}", file=sys.stderr)

    if not processes:
        print("No valid status files found", file=sys.stderr)
        sys.exit(1)

    if not submitter_analysis_id:
        print("No submitter_analysis_id found in status files", file=sys.stderr)
        sys.exit(1)

    # Sort processes by timestamp
    processes.sort(key=lambda x: x['timestamp'])

    # Generate individual receipt
    generate_individual_receipt(processes, analysis_info, args.output_json)

    print(f"Generated individual receipt for analysis {submitter_analysis_id} with {len(processes)} processes", file=sys.stderr)
    print(f"JSON: {args.output_json}", file=sys.stderr)


if __name__ == '__main__':
    main()
