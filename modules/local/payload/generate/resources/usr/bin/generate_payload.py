#!/usr/bin/env python3

import json
import sys
import argparse
import csv
from pathlib import Path

def detect_delimiter(file_path):
    """Detect delimiter by examining the first few lines of the file"""
    try:
        with open(file_path, 'r', newline='', encoding='utf-8') as f:
            # Read first few lines to get a better sample
            lines = []
            for i, line in enumerate(f):
                if i >= 3:  # Check first 3 lines
                    break
                lines.append(line.strip())
            
            if not lines:
                return '\t'  # Default to tab
            
            # Count occurrences of common delimiters across all lines
            tab_count = sum(line.count('\t') for line in lines)
            comma_count = sum(line.count(',') for line in lines)
            
            # # Additional check: if file extension suggests format
            # file_ext = Path(file_path).suffix.lower()
            # if file_ext == '.csv':
            #     return ','
            # elif file_ext == '.tsv':
            #     return '\t'
            
            # If there are significantly more commas than tabs, likely CSV
            # If there are more tabs than commas, likely TSV
            # Use a small threshold to handle cases where there might be commas in quoted fields
            if comma_count > tab_count and comma_count > 0:
                return ','
            else:
                return '\t'
                
    except Exception:
        return '\t'  # Default to tab if detection fails

def read_metadata_file(file_path):
    """Read metadata file (TSV or CSV) and return list of dictionaries with headers as keys"""
    try:
        delimiter = detect_delimiter(file_path)
        file_type = "CSV" if delimiter == ',' else "TSV"
        
        with open(file_path, 'r', newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            data = list(reader)
            
        print(f"Successfully read {file_type} file {file_path} with {len(data)} rows")
        return data
        
    except Exception as e:
        print(f"Error reading metadata file {file_path}: {e}", file=sys.stderr)
        return []

def main():
    parser = argparse.ArgumentParser(description='Generate JSON payload from metadata files (TSV or CSV format)')
    parser.add_argument('--submitter-analysis-id', required=True, help='Submitter analysis ID')
    parser.add_argument('--analysis-type', required=True, help='Analysis type')
    parser.add_argument('--study-id', required=True, help='Study ID')
    parser.add_argument('--file-meta', required=True, help='Path to file metadata (TSV or CSV)')
    parser.add_argument('--analysis-meta', required=True, help='Path to analysis metadata (TSV or CSV)')
    parser.add_argument('--workflow-meta', help='Path to workflow metadata (TSV or CSV)')
    parser.add_argument('--data-files', nargs='*', help='Paths to data files')
    parser.add_argument('--output', '-o', required=True, help='Output JSON file path')
    
    args = parser.parse_args()
    
    # Validate required metadata based on analysis type
    if args.analysis_type != "sequenceExperiment" and not args.workflow_meta:
        print(f"Error: workflow-meta is required for analysis type '{args.analysis_type}'", file=sys.stderr)
        sys.exit(1)
    
    # Read metadata files (TSV or CSV format)
    file_data = read_metadata_file(args.file_meta)
    if not file_data:
        print(f"Error: Could not read file metadata from {args.file_meta}", file=sys.stderr)
        sys.exit(1)
    
    analysis_data = read_metadata_file(args.analysis_meta)
    if not analysis_data:
        print(f"Error: Could not read analysis metadata from {args.analysis_meta}", file=sys.stderr)
        sys.exit(1)
    
    # Convert list of rows to single dict for analysis_meta (assuming single row)
    analysis_metadata = analysis_data[0] if analysis_data else {}
    
    workflow_metadata = None
    if args.workflow_meta:
        workflow_data = read_metadata_file(args.workflow_meta)
        if not workflow_data:
            print(f"Error: Could not read workflow metadata from {args.workflow_meta}", file=sys.stderr)
            sys.exit(1)
        # Convert list of rows to single dict for workflow_meta (assuming single row)
        workflow_metadata = workflow_data[0] if workflow_data else {}
    
    # Process files from file metadata (TSV or CSV format)
    files_info = []
    for file_row in file_data:
        file_info = {
            "fileName": file_row.get("fileName", ""),
            "fileSize": int(file_row.get("fileSize", 0)) if file_row.get("fileSize", "").isdigit() else 0,
            "dataType": file_row.get("dataType", ""),
            "fileAccess": file_row.get("fileAccess", "controlled"),
            "fileMd5sum": file_row.get("fileMd5sum", ""),
            "fileType": file_row.get("fileType", "")
        }
        files_info.append(file_info)
    
    # Create payload structure based on analysis type
    payload = {
        "analysisType": {
            "name": args.analysis_type
        },
        "studyId": args.study_id,
        "submitter_analysis_id": args.submitter_analysis_id,
        "files": files_info
    }
    
    # Add analysis metadata fields from TSV data
    if "data_category" in analysis_metadata:
        payload["data_category"] = analysis_metadata["data_category"]
    if "submitter_experiment_id" in analysis_metadata:
        payload["submitter_experiment_id"] = analysis_metadata["submitter_experiment_id"]
    
    # Add fields specific to analysis type
    if args.analysis_type == "sequenceExperiment":
        # For sequenceExperiment, only file_meta and analysis_meta are used
        pass  # Base payload is sufficient
        
    elif args.analysis_type == "sequenceAlignment":
        # Add genome_build and workflow for sequenceAlignment
        if "genome_build" in analysis_metadata:
            payload["genome_build"] = analysis_metadata["genome_build"]
        
        if workflow_metadata:
            payload["workflow"] = {
                "submitter_workflow_id": workflow_metadata.get("submitter_workflow_id", ""),
                "workflow_name": workflow_metadata.get("workflow_name", ""),
                "workflow_url": workflow_metadata.get("workflow_url", ""),
                "workflow_version": workflow_metadata.get("workflow_version", "")
            }
            
    elif args.analysis_type == "variantCall":
        # Add genome_annotation, genome_build, variant fields and workflow for variantCall
        if "genome_annotation" in analysis_metadata:
            payload["genome_annotation"] = analysis_metadata["genome_annotation"]
        if "genome_build" in analysis_metadata:
            payload["genome_build"] = analysis_metadata["genome_build"]
        if "variant_calling_strategy" in analysis_metadata:
            payload["variant_calling_strategy"] = analysis_metadata["variant_calling_strategy"]
        if "variant_class" in analysis_metadata:
            payload["variant_class"] = analysis_metadata["variant_class"]
        
        if workflow_metadata:
            payload["workflow"] = {
                "submitter_workflow_id": workflow_metadata.get("submitter_workflow_id", ""),
                "workflow_name": workflow_metadata.get("workflow_name", ""),
                "workflow_url": workflow_metadata.get("workflow_url", ""),
                "workflow_version": workflow_metadata.get("workflow_version", "")
            }
    
    # Write payload to JSON file
    with open(args.output, "w") as f:
        json.dump(payload, f, indent=2)
    
    print(f"Generated payload: {args.output}")

if __name__ == "__main__":
    main()
