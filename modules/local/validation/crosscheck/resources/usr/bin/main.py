#!/usr/bin/env python3
"""
MD5 checksum validation script for genomic data submission workflow.
Cross-validates actual file MD5 checksums against payload metadata.
"""

import json
import sys
import os
import hashlib
import argparse
from pathlib import Path


def calculate_md5(file_path):
    """Calculate MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    try:
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    except Exception as e:
        print(f'ERROR: Failed to calculate MD5 for {file_path}: {e}', file=sys.stderr)
        return None


def validate_md5_checksums(payload_file, files_list):
    """Validate MD5 checksums of files against payload metadata."""
    try:
        # Load analysis payload
        with open(payload_file, 'r') as f:
            payload = json.load(f)
        
        print(f'Loaded payload for MD5 validation: {payload.get("studyId", "unknown")}')
        print(f'Files to validate: {len(files_list)} files')
        
        # Create payload file mapping by filename
        payload_files = payload.get('files', [])
        file_md5_map = {}
        
        for pf in payload_files:
            if 'fileName' in pf and 'fileMd5sum' in pf:
                file_md5_map[pf['fileName']] = pf['fileMd5sum']
        
        print(f'Found {len(file_md5_map)} files with MD5 checksums in payload')
        
        # Validate MD5 checksums for each file
        validation_errors = []
        files_validated = 0
        
        for actual_file in files_list:
            filename = os.path.basename(actual_file)
            
            if filename not in file_md5_map:
                validation_errors.append(f'File {filename} not found in payload metadata')
                continue
            
            expected_md5 = file_md5_map[filename]
            if not expected_md5:
                validation_errors.append(f'File {filename} has no MD5 checksum in payload')
                continue
            
            print(f'Validating MD5 for {filename}...')
            actual_md5 = calculate_md5(actual_file)
            
            if actual_md5 is None:
                validation_errors.append(f'Failed to calculate MD5 for {filename}')
                continue
            
            if actual_md5.lower() != expected_md5.lower():
                validation_errors.append(f'MD5 mismatch for {filename}: expected {expected_md5}, got {actual_md5}')
            else:
                print(f'MD5 validation passed for {filename}')
                files_validated += 1
        
        if validation_errors:
            print('MD5 VALIDATION ERRORS:', file=sys.stderr)
            for error in validation_errors:
                print(f'  - {error}', file=sys.stderr)
            return False, validation_errors
        else:
            print(f'MD5 checksum validation completed successfully for {files_validated} files')
            return True, []

    except json.JSONDecodeError as e:
        error_msg = f'Invalid JSON in payload file: {e}'
        print(f'ERROR: {error_msg}', file=sys.stderr)
        return False, [error_msg]
    except Exception as e:
        error_msg = f'MD5 checksum validation failed: {e}'
        print(f'ERROR: {error_msg}', file=sys.stderr)
        return False, [error_msg]


def main():
    """Main function for command line execution."""
    parser = argparse.ArgumentParser(description='Validate MD5 checksums of files against analysis payload metadata')
    parser.add_argument('payload_file', help='Path to analysis payload JSON file')
    parser.add_argument('files', nargs='+', help='List of files to validate')
    
    args = parser.parse_args()
    
    try:
        success, errors = validate_md5_checksums(args.payload_file, args.files)
        
        if success:
            sys.exit(0)
        else:
            sys.exit(1)
    
    except Exception as e:
        print(f'Unexpected error during MD5 validation: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
