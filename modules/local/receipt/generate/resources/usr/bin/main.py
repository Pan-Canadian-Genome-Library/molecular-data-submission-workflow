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


def parse_analysis_file(file_path):
    """Parse analysis JSON file and extract required information."""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        return {
            'song_analysis_id': data.get('analysisId', 'unknown'),
            'submitter_analysis_id': data.get('submitter_analysis_id', 'unknown'),
            'analysis_state': data.get('analysisState', 'unknown'),
            'published_at': data.get('publishedAt', 'unknown'),
            'study_id': data.get('studyId', 'unknown'),
            'analysis_type': data.get('analysisType', {}).get('name', 'unknown') if isinstance(data.get('analysisType'), dict) else 'unknown'
        }
    except Exception as e:
        print(f"Error parsing analysis file {file_path}: {e}", file=sys.stderr)
        return {
            'song_analysis_id': 'unknown',
            'submitter_analysis_id': 'unknown',
            'analysis_state': 'unknown', 
            'published_at': 'unknown',
            'study_id': 'unknown',
            'analysis_type': 'unknown'
        }


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
        
        # Extract error message from various possible locations
        error_message = ''
        if 'error_message' in data:
            error_message = data['error_message']
        elif 'details' in data and 'error_message' in data['details']:
            error_message = data['details']['error_message']
        elif 'details' in data and 'error_details' in data['details']:
            error_message = data['details']['error_details']
        
        return {
            'submitter_analysis_id': submitter_analysis_id,
            'process': data.get('process', 'UNKNOWN'),
            'status': data.get('status', 'UNKNOWN'),
            'exit_code': data.get('exit_code', -1),
            'timestamp': data.get('timestamp', ''),
            'details': data.get('details', {}),
            'error_message': error_message.strip() if error_message else '',
            'source_file': str(file_path)
        }
    except Exception as e:
        print(f"Error parsing {file_path}: {e}", file=sys.stderr)
        return None


def generate_individual_receipt(processes, analysis_info, output_file):
    """Generate individual analysis receipt in JSON format."""
    # Determine overall status for this analysis
    overall_status = 'SUCCESS'
    if any(p['status'] == 'FAILED' for p in processes):
        overall_status = 'FAILED'
    
    # Clean up process data for JSON output
    clean_processes = []
    for p in processes:
        clean_process = {
            'process': p['process'],
            'status': p['status'],
            'exit_code': p['exit_code'],
            'timestamp': p['timestamp']
        }
        if p['error_message']:
            clean_process['error_message'] = p['error_message']
        if p['details']:
            clean_process['details'] = p['details']
        clean_processes.append(clean_process)
    
    # Create individual receipt structure
    receipt = {
        'submitter_analysis_id': analysis_info.get('submitter_analysis_id', processes[0]['submitter_analysis_id']),
        'song_analysis_id': analysis_info.get('song_analysis_id', 'unknown'),
        'overall_status': overall_status,
        'analysis_type': analysis_info.get('analysis_type', 'unknown'),
        'study_id': analysis_info.get('study_id', 'unknown'),
        'analysis_state': analysis_info.get('song_analysis_state', 'unknown'),
        'published_at': analysis_info.get('song_analysis_publish_at', 'unknown'),
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
                        help='Analysis JSON file (optional, only for song_analysis_id, song_analysis_state, song_analysis_publish_at)')
    parser.add_argument('--submitter-analysis-id', required=True, help='Submitter analysis id (from meta)')
    parser.add_argument('--study-id', required=True, help='Study id (from meta)')
    parser.add_argument('--analysis-type', required=True, help='Analysis type (from meta)')
    parser.add_argument('--output-json', required=True,
                        help='Output JSON file path')

    args = parser.parse_args()

    # Parse analysis file if provided, else use defaults for song fields
    song_analysis_id = 'unknown'
    song_analysis_state = 'unknown'
    song_analysis_publish_at = 'unknown'
    if args.analysis_file:
        try:
            with open(args.analysis_file, 'r') as f:
                data = json.load(f)
            song_analysis_id = data.get('analysisId', 'unknown')
            song_analysis_state = data.get('analysisState', 'unknown')
            song_analysis_publish_at = data.get('publishedAt', 'unknown')
        except Exception as e:
            print(f"Error parsing analysis file {args.analysis_file}: {e}", file=sys.stderr)
    else:
        print("No analysis file provided, using default values for song_analysis_id, song_analysis_state, song_analysis_publish_at.", file=sys.stderr)

    # Build analysis_info dict for receipt
    analysis_info = {
        'song_analysis_id': song_analysis_id,
        'song_analysis_state': song_analysis_state,
        'song_analysis_publish_at': song_analysis_publish_at,
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
        if data and data['submitter_analysis_id']:
            if submitter_analysis_id is None:
                submitter_analysis_id = data['submitter_analysis_id']
            elif submitter_analysis_id != data['submitter_analysis_id']:
                print(f"Warning: Mixed submitter_analysis_ids detected: {submitter_analysis_id} vs {data['submitter_analysis_id']}", file=sys.stderr)
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
