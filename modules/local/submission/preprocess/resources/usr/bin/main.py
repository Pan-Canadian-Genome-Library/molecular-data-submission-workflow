#!/usr/bin/env python3

import csv
import json
import os
import sys
import yaml
from pathlib import Path

def filter_by_analysis_id(filepath, id_column, target_id):
    """Filter CSV/TSV file by analysis ID and return matching rows"""
    filtered_rows = []
    try:
        # Try TSV first, then CSV
        for delimiter in ['\t', ',']:
            try:
                with open(filepath, 'r') as f:
                    reader = csv.DictReader(f, delimiter=delimiter)
                    for row in reader:
                        if row.get(id_column) == target_id:
                            filtered_rows.append(row)
                if filtered_rows:  # If we found rows, we used the right delimiter
                    break
            except:
                continue
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return filtered_rows

def main():
    # Get command line arguments
    if len(sys.argv) != 13:
        print("Usage: preprocess_submission.py <analysis_id> <study_id> <analysis_type> <status> <file_metadata> <analysis_metadata> <workflow_metadata> <read_group_metadata> <experiment_metadata> <specimen_metadata> <sample_metadata> <files_directory>")
        sys.exit(1)
    
    analysis_id = sys.argv[1]
    study_id = sys.argv[2]
    analysis_type = sys.argv[3]
    status = sys.argv[4]
    file_metadata_path = sys.argv[5]
    analysis_metadata_path = sys.argv[6]
    workflow_metadata_path = sys.argv[7]
    read_group_metadata_path = sys.argv[8]
    experiment_metadata_path = sys.argv[9]
    specimen_metadata_path = sys.argv[10]
    sample_metadata_path = sys.argv[11]
    files_directory = sys.argv[12]
    
    # Check if we have real input files
    use_real_data = file_metadata_path and file_metadata_path != "NO_FILE" and os.path.exists(file_metadata_path)
    
    if not use_real_data:
        print(f"Error: Input file {file_metadata_path} does not exist or is invalid")
        sys.exit(1)
    
    # Extract file metadata for this analysis
    file_rows = filter_by_analysis_id(file_metadata_path, 'submitter_analysis_id', analysis_id)
    with open(f'file_meta_{analysis_id}.tsv', 'w', newline='') as f:
        if file_rows:
            writer = csv.DictWriter(f, fieldnames=file_rows[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(file_rows)
        else:
            # Fallback to empty file with headers
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['submitter_analysis_id', 'fileName', 'fileSize', 'fileMd5sum', 'fileType', 'fileAccess', 'dataType'])
    
    # Extract analysis metadata
    analysis_rows = filter_by_analysis_id(analysis_metadata_path, 'submitter_analysis_id', analysis_id)
    with open(f'analysis_meta_{analysis_id}.tsv', 'w', newline='') as f:
        if analysis_rows:
            writer = csv.DictWriter(f, fieldnames=analysis_rows[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(analysis_rows)
        else:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['studyId', 'submitter_analysis_id', 'analysisType', 'submitter_participant_id', 'submitter_specimen_id', 'submitter_sample_id', 'submitter_experiment_id'])
    
    # Extract workflow metadata
    workflow_rows = filter_by_analysis_id(workflow_metadata_path, 'submitter_analysis_id', analysis_id)
    with open(f'workflow_meta_{analysis_id}.tsv', 'w', newline='') as f:
        if workflow_rows:
            writer = csv.DictWriter(f, fieldnames=workflow_rows[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(workflow_rows)
        else:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['submitter_workflow_id', 'submitter_analysis_id', 'workflow_name', 'workflow_version', 'workflow_url'])
    
    # Process other metadata files similarly
    for metadata_file, output_file, id_col in [
        (specimen_metadata_path, f'specimen_{analysis_id}.tsv', 'submitter_specimen_id'),
        (sample_metadata_path, f'sample_{analysis_id}.tsv', 'submitter_sample_id'),
        (experiment_metadata_path, f'experiment_{analysis_id}.tsv', 'submitter_experiment_id'),
        (read_group_metadata_path, f'read_group_{analysis_id}.tsv', 'submitter_read_group_id')
    ]:
        # For these files, we need to cross-reference with analysis metadata to find the right IDs
        entity_rows = []
        if analysis_rows:
            for analysis_row in analysis_rows:
                if id_col == 'submitter_specimen_id':
                    target_id = analysis_row.get('submitter_specimen_id')
                elif id_col == 'submitter_sample_id':
                    target_id = analysis_row.get('submitter_sample_id')
                elif id_col == 'submitter_experiment_id':
                    target_id = analysis_row.get('submitter_experiment_id')
                else:  # read_group
                    target_id = analysis_row.get('submitter_experiment_id')  # Read groups are linked to experiments
                    id_col = 'submitter_experiment_id'
                    
                if target_id:
                    entity_rows.extend(filter_by_analysis_id(metadata_file, id_col, target_id))
        
        with open(output_file, 'w', newline='') as f:
            if entity_rows:
                writer = csv.DictWriter(f, fieldnames=entity_rows[0].keys(), delimiter='\t')
                writer.writeheader()
                writer.writerows(entity_rows)
            else:
                # Create empty file with basic headers
                writer = csv.writer(f, delimiter='\t')
                if 'specimen' in output_file:
                    writer.writerow(['submitter_participant_id', 'submitter_specimen_id'])
                elif 'sample' in output_file:
                    writer.writerow(['submitter_sample_id', 'submitter_specimen_id'])
                elif 'experiment' in output_file:
                    writer.writerow(['submitter_experiment_id', 'submitter_sample_id'])
                elif 'read_group' in output_file:
                    writer.writerow(['submitter_read_group_id', 'submitter_experiment_id'])
    
    # Create data files list from files directory and file metadata
    data_files_list = []
    if file_rows and files_directory and files_directory != "NO_DIR":
        for file_row in file_rows:
            filename = file_row.get('fileName', '')
            if filename:
                # Check if files_directory is a staged directory or path
                if os.path.isdir(files_directory):
                    # If it's a directory, look for the file inside it
                    file_path = os.path.join(files_directory, filename)
                else:
                    # If it's not a directory, assume it's a path string
                    file_path = os.path.join(str(files_directory), filename)
                
                if os.path.exists(file_path):
                    data_files_list.append(os.path.abspath(file_path))
                else:
                    # Also try looking in the current working directory
                    # in case files were staged directly
                    alt_path = filename
                    if os.path.exists(alt_path):
                        data_files_list.append(os.path.abspath(alt_path))
    
    with open(f'data_files_{analysis_id}.txt', 'w') as f:
        for file_path in data_files_list:
            f.write(f"{file_path}\n")

if __name__ == "__main__":
    main()
