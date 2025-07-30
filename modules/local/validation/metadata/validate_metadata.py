#!/usr/bin/env python3
"""
Metadata validation script for genomic data submission workflow.
Validates consistency between analysis payload and read group data.
"""

import json
import sys
import csv
import os
import argparse
from pathlib import Path


def load_analysis_payload(payload_file):
    """Load and validate basic structure of analysis payload JSON."""
    try:
        with open(payload_file, 'r') as f:
            payload = json.load(f)
        
        # Check required fields
        required_fields = ['studyId', 'analysisType', 'files']
        missing_fields = []
        
        for field in required_fields:
            if field not in payload:
                missing_fields.append(field)
        
        if missing_fields:
            raise ValueError(f'Missing required fields: {", ".join(missing_fields)}')
        
        return payload
    
    except FileNotFoundError:
        raise FileNotFoundError(f'Analysis payload file not found: {payload_file}')
    except json.JSONDecodeError as e:
        raise ValueError(f'Invalid JSON in payload file: {e}')


def load_read_group_data(read_group_file):
    """Load read group data from TSV file."""
    try:
        read_groups = []
        with open(read_group_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            read_groups = list(reader)
        
        if not read_groups:
            raise ValueError('No read groups found in read group file')
        
        return read_groups
    
    except FileNotFoundError:
        raise FileNotFoundError(f'Read group file not found: {read_group_file}')
    except Exception as e:
        raise ValueError(f'Error reading read group file: {e}')


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
        rg_id = rg.get('submitter_read_group_id', 'unknown')
        
        # Check file_r1
        file_r1 = rg.get('file_r1', '').strip()
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
        library_layout = rg.get('library_layout', '').upper()
        if library_layout == 'PAIRED':
            if not file_r1:
                validation_errors.append(f'Read group {rg_id}: library_layout is PAIRED but file_r1 is missing')
            if not file_r2:
                validation_errors.append(f'Read group {rg_id}: library_layout is PAIRED but file_r2 is missing')
        elif library_layout == 'SINGLE':
            if not file_r1:
                validation_errors.append(f'Read group {rg_id}: library_layout is SINGLE but file_r1 is missing')
            if file_r2:
                validation_errors.append(f'Read group {rg_id}: library_layout is SINGLE but file_r2 is provided')
        elif library_layout and library_layout not in ['PAIRED', 'SINGLE']:
            validation_errors.append(f'Read group {rg_id}: invalid library_layout "{library_layout}", must be PAIRED or SINGLE')
    
    print(f'Found {len(read_group_files)} unique files referenced in read groups: {sorted(read_group_files)}')
    
    # Check for payload files not referenced in read groups
    unreferenced_files = payload_filenames - read_group_files
    if unreferenced_files:
        for uf in sorted(unreferenced_files):
            validation_errors.append(f'Payload file "{uf}" not referenced in any read group')
    
    return validation_errors


def validate_read_group_fields(read_groups):
    """Validate required fields in read groups."""
    validation_errors = []
    
    # Check for additional required read group fields
    required_rg_fields = ['submitter_read_group_id', 'submitter_experiment_id', 'library_name', 'library_layout']
    
    for i, rg in enumerate(read_groups):
        rg_id = rg.get('submitter_read_group_id', f'read_group_{i+1}')
        
        for field in required_rg_fields:
            if not rg.get(field, '').strip():
                validation_errors.append(f'Read group {rg_id}: missing required field "{field}"')
        
        # Additional validation for specific fields
        platform_unit = rg.get('platform_unit', '').strip()
        if platform_unit and len(platform_unit) > 100:
            validation_errors.append(f'Read group {rg_id}: platform_unit too long (max 100 characters)')
        
        # Validate read lengths if provided
        try:
            read_length_r1 = rg.get('read_length_r1', '').strip()
            if read_length_r1:
                length_r1 = int(read_length_r1)
                if length_r1 <= 0:
                    validation_errors.append(f'Read group {rg_id}: read_length_r1 must be positive integer')
        except ValueError:
            validation_errors.append(f'Read group {rg_id}: read_length_r1 must be a valid integer')
        
        try:
            read_length_r2 = rg.get('read_length_r2', '').strip()
            if read_length_r2:
                length_r2 = int(read_length_r2)
                if length_r2 <= 0:
                    validation_errors.append(f'Read group {rg_id}: read_length_r2 must be positive integer')
        except ValueError:
            validation_errors.append(f'Read group {rg_id}: read_length_r2 must be a valid integer')
        
        # Validate insert size if provided
        try:
            insert_size = rg.get('insert_size', '').strip()
            if insert_size:
                size = int(insert_size)
                if size <= 0:
                    validation_errors.append(f'Read group {rg_id}: insert_size must be positive integer')
        except ValueError:
            validation_errors.append(f'Read group {rg_id}: insert_size must be a valid integer')
    
    return validation_errors


def validate_payload_analysis_type(payload):
    """Validate analysis type structure in payload."""
    validation_errors = []
    
    analysis_type = payload.get('analysisType', {})
    if not isinstance(analysis_type, dict):
        validation_errors.append('analysisType must be an object')
        return validation_errors
    
    # Check for required analysis type fields
    required_analysis_fields = ['name', 'version']
    for field in required_analysis_fields:
        if field not in analysis_type:
            validation_errors.append(f'Missing required analysis type field: {field}')
    
    return validation_errors


def basic_payload_validation(payload_file):
    """Perform basic payload validation without read group data."""
    try:
        payload = load_analysis_payload(payload_file)
        validation_errors = validate_payload_analysis_type(payload)
        
        if validation_errors:
            return False, validation_errors
        
        print('Basic metadata validation completed successfully')
        return True, []
    
    except Exception as e:
        return False, [str(e)]


def comprehensive_validation(payload_file, read_group_file):
    """Perform comprehensive validation with both payload and read group data."""
    try:
        # Load data
        payload = load_analysis_payload(payload_file)
        read_groups = load_read_group_data(read_group_file)
        
        print(f'Loaded {len(read_groups)} read groups from biospecimen data')
        
        # Perform all validations
        validation_errors = []
        validation_errors.extend(validate_payload_analysis_type(payload))
        validation_errors.extend(validate_read_group_fields(read_groups))
        validation_errors.extend(validate_read_group_files(payload, read_groups))
        
        if validation_errors:
            return False, validation_errors
        
        print('Read group file validation completed successfully')
        print('All file references are consistent between read groups and payload')
        return True, []
    
    except Exception as e:
        return False, [str(e)]


def main():
    """Main function for command line execution."""
    parser = argparse.ArgumentParser(description='Validate metadata consistency between payload and read group data')
    parser.add_argument('payload_file', help='Path to analysis payload JSON file')
    parser.add_argument('--read-group-file', help='Path to read group TSV file (optional)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Set verbosity
    if not args.verbose:
        # Redirect stdout to stderr for non-verbose mode to keep validation output
        pass
    
    try:
        if args.read_group_file and os.path.exists(args.read_group_file):
            success, errors = comprehensive_validation(args.payload_file, args.read_group_file)
        else:
            if args.read_group_file:
                print(f"Warning: Read group file not found: {args.read_group_file}")
                print("Falling back to basic payload validation")
            success, errors = basic_payload_validation(args.payload_file)
        
        if not success:
            print('VALIDATION ERRORS:')
            for error in errors:
                print(f'  - {error}')
            sys.exit(1)
        else:
            print('Metadata validation completed successfully')
            sys.exit(0)
    
    except Exception as e:
        print(f'Unexpected error during validation: {e}')
        sys.exit(1)


if __name__ == '__main__':
    main()
