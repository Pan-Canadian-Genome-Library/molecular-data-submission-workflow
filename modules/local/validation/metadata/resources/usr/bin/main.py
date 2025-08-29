#!/usr/bin/env python3
"""
Metadata cross-validation script for genomic data submission workflow.
Cross-validates file consistency between analysis payload and read group data.
Assumes upstream validation of individual files and biospecimen relationships.
"""

import json
import sys
import csv
import os
import argparse
from pathlib import Path


def load_analysis_payload(payload_file):
    """Load analysis payload JSON (assumes upstream validation)."""
    try:
        with open(payload_file, 'r') as f:
            payload = json.load(f)
        
        print(f"Loaded analysis payload with {len(payload.get('files', []))} files")
        return payload
    
    except Exception as e:
        print(f'Error loading payload file: {e}', file=sys.stderr)
        raise


def load_biospecimen_data(file_path, file_type):
    """Load biospecimen data from TSV file (assumes upstream validation)."""
    try:
        data = []
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            data = list(reader)
        
        print(f"Loaded {len(data)} {file_type} records")
        return data
    
    except Exception as e:
        print(f"Error loading {file_type} file {file_path}: {e}", file=sys.stderr)
        raise


def validate_read_group_files(payload, read_groups):
    """Validate file references between payload and read groups."""
    validation_errors = []
    
    # Get payload files
    payload_files = payload.get('files', [])
    payload_filenames = set()
    
    for file_entry in payload_files:
        if 'fileName' in file_entry:
            payload_filenames.add(file_entry['fileName'])
    
    print(f'Found {len(payload_filenames)} files in payload: {sorted(payload_filenames)}')
    
    # Cross-validate read group file references with payload files
    read_group_files = set()
    
    for rg in read_groups:
        rg_id = rg.get('submitter_read_group_id')
        
        # Check file_r1
        file_r1 = rg.get('file_r1').strip()
        if file_r1:
            read_group_files.add(file_r1)
            if file_r1 not in payload_filenames:
                validation_errors.append(f'Read group {rg_id}: file_r1 "{file_r1}" not found in payload files')
        
        # Check file_r2
        file_r2 = rg.get('file_r2', '').strip()
        if file_r2:
            read_group_files.add(file_r2)
            if file_r2 not in payload_filenames:
                validation_errors.append(f'Read group {rg_id}: file_r2 "{file_r2}" not found in payload files')
        
        # Validate library layout consistency
        library_layout = rg.get('library_layout').upper()
        if 'PAIRED' in library_layout:
            if not file_r1:
                validation_errors.append(f'Read group {rg_id}: library_layout is PAIRED but file_r1 is missing')
            if not file_r2:
                validation_errors.append(f'Read group {rg_id}: library_layout is PAIRED but file_r2 is missing')
        elif 'SINGLE' in library_layout:
            if not file_r1:
                validation_errors.append(f'Read group {rg_id}: library_layout is SINGLE but file_r1 is missing')
            if file_r2:
                validation_errors.append(f'Read group {rg_id}: library_layout is SINGLE but file_r2 is provided')
        else:
            pass
    

    print(f'Found {len(read_group_files)} files referenced in read groups: {sorted(read_group_files)}')
  
    return validation_errors


def comprehensive_validation(payload_file, specimen_file, sample_file, experiment_file, read_group_file, analysis_type):
    """Perform cross-validation between payload and biospecimen entity data (assumes upstream validation)."""
    try:
        # Load and validate payload
        payload = load_analysis_payload(payload_file)
        
        # Extract analysis type from payload if not provided
        if not analysis_type and 'analysisType' in payload:
            analysis_type = payload['analysisType'].get('name', '')
        
        print(f'Analysis type: {analysis_type}')
        
        # Check if read group is required (only for sequenceExperiment)
        read_group_required = analysis_type == "sequenceExperiment"
        
        # Disable the loading of biospecimen entities
        # # Load biospecimen data
        # specimens = load_biospecimen_data(specimen_file, 'specimen')
        # samples = load_biospecimen_data(sample_file, 'sample')
        # experiments = load_biospecimen_data(experiment_file, 'experiment')
        
        # Load read groups if required or if file is provided
        read_groups = []
        if read_group_required or (read_group_file and os.path.exists(read_group_file) and os.path.getsize(read_group_file) > 0):
            read_groups = load_biospecimen_data(read_group_file, 'read_group')
        
        # Validate read group requirement
        if read_group_required and not read_groups:
            raise ValueError(f'Read group data is required for analysis type "{analysis_type}" but no valid read group data found')
        
        print(f'Read group required: {read_group_required}')
        # print(f'Loaded biospecimen data: {len(specimens)} specimens, {len(samples)} samples, {len(experiments)} experiments, {len(read_groups)} read groups')
        
        # Perform cross-validation checks (individual validation done upstream)
        validation_errors = []
        
        # Cross-validate read group file references with payload (if read groups exist)
        if read_groups:
            validation_errors.extend(validate_read_group_files(payload, read_groups))
        
        if validation_errors:
            return False, validation_errors
        
        return True, []
    
    except Exception as e:
        # All exceptions from comprehensive_validation are sent to stderr
        return False, [f'Unexpected validation error: {str(e)}']


def main():
    """Main function for command line execution."""
    parser = argparse.ArgumentParser(description='Cross-validate file consistency between payload and read group data (assumes upstream validation of individual files and relationships)')
    parser.add_argument('payload_file', help='Path to analysis payload JSON file')
    parser.add_argument('--specimen-file', required=False, help='Path to specimen TSV file')
    parser.add_argument('--sample-file', required=False, help='Path to sample TSV file')
    parser.add_argument('--experiment-file', required=False, help='Path to experiment TSV file')
    parser.add_argument('--read-group-file', required=False, help='Path to read group TSV file (analysis type dependent)')
    parser.add_argument('--analysis-type', required=True, help='Analysis type from metadata (optional, will be extracted from payload if available)')
    
    args = parser.parse_args()
    
    try:
        # Run comprehensive validation with all data
        success, errors = comprehensive_validation(
            args.payload_file, 
            args.specimen_file, 
            args.sample_file, 
            args.experiment_file, 
            args.read_group_file,
            args.analysis_type
        )
        
        if not success:
            # Send error messages to stderr with consistent format
            print('VALIDATION ERRORS:', file=sys.stderr)
            for error in errors:
                print(f'  - {error}', file=sys.stderr)
            sys.exit(1)
        else:
            print('Payload-Biospecimen files consistency validation completed successfully')
            sys.exit(0)
    
    except Exception as e:
        print(f'Unexpected error during validation: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
