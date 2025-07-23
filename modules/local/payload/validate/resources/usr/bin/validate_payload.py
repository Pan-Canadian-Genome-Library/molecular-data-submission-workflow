#!/usr/bin/env python3

import json
import sys
import argparse
import urllib.request
import urllib.error
from pathlib import Path

def download_schema(schema_url):
    """Download JSON schema from URL"""
    try:
        with urllib.request.urlopen(schema_url) as response:
            schema_content = response.read().decode('utf-8')
            return json.loads(schema_content)
    except urllib.error.URLError as e:
        print(f"Error downloading schema from {schema_url}: {e}", file=sys.stderr)
        return None
    except json.JSONDecodeError as e:
        print(f"Error parsing schema JSON: {e}", file=sys.stderr)
        return None

def load_payload(payload_file):
    """Load JSON payload from file"""
    try:
        with open(payload_file, 'r', encoding='utf-8') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error parsing payload JSON: {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading payload file: {e}", file=sys.stderr)
        return None

def validate_payload(payload, schema):
    """Validate payload against schema using jsonschema library"""
    try:
        import jsonschema
    except ImportError:
        # Try common package locations first
        import subprocess
        
        # Check if jsonschema is available under different import paths
        for module_name in ['jsonschema', 'jsonschema.validate']:
            try:
                __import__(module_name)
                import jsonschema
                break
            except ImportError:
                continue
        else:
            # If still not found, provide a clear error message
            error_msg = ("jsonschema library not available in container. "
                        "Please ensure the container includes jsonschema or update the container image.")
            print(error_msg, file=sys.stderr)
            return False, error_msg
    
    # Extract the actual schema from the wrapper if it exists
    # The downloaded schema may be wrapped in metadata, we need the actual JSON schema
    actual_schema = schema
    if isinstance(schema, dict) and 'schema' in schema:
        actual_schema = schema['schema']
        print(f"Extracted actual schema from wrapper", file=sys.stderr)
    else:
        print(f"Using schema as-is (no wrapper detected)", file=sys.stderr)
    
    try:
        jsonschema.validate(instance=payload, schema=actual_schema)
        return True, "Validation successful"
    except jsonschema.ValidationError as e:
        return False, f"Validation failed: {e.message}"
    except jsonschema.SchemaError as e:
        return False, f"Schema error: {e.message}"
    except Exception as e:
        return False, f"Validation error: {str(e)}"

def main():
    parser = argparse.ArgumentParser(description='Validate JSON payload against schema from URL')
    parser.add_argument('--payload', required=True, help='Path to JSON payload file')
    parser.add_argument('--schema-url', required=True, help='URL to download JSON schema')
    parser.add_argument('--output-report', required=True, help='Path to write validation report')
    
    args = parser.parse_args()
    
    # Load payload first to get analysisType
    payload = load_payload(args.payload)
    if payload is None:
        with open(args.output_report, 'w') as f:
            f.write("STATUS: INVALID\n")
            f.write("MESSAGE: Failed to load payload file\n")
            f.write(f"PAYLOAD_FILE: {args.payload}\n")
        print("Failed to load payload file", file=sys.stderr)
        sys.exit(1)  # Script execution error, not validation error
    
    # Extract analysisType from payload to construct final schema URL
    analysis_type = payload.get('analysisType', {}).get('name', '')
    if not analysis_type:
        with open(args.output_report, 'w') as f:
            f.write("STATUS: INVALID\n")
            f.write("MESSAGE: analysisType.name not found in payload\n")
            f.write(f"PAYLOAD_FILE: {args.payload}\n")
        print("analysisType.name not found in payload", file=sys.stderr)
        sys.exit(1)  # Script execution error, not validation error
    
    # Construct final schema URL
    final_schema_url = f"{args.schema_url}/{analysis_type}"
    
    # Download schema
    schema = download_schema(final_schema_url)
    if schema is None:
        with open(args.output_report, 'w') as f:
            f.write("STATUS: INVALID\n")
            f.write(f"MESSAGE: Failed to download schema from {final_schema_url}\n")
            f.write(f"PAYLOAD_FILE: {args.payload}\n")
            f.write(f"ANALYSIS_TYPE: {analysis_type}\n")
            f.write(f"BASE_SCHEMA_URL: {args.schema_url}\n")
            f.write(f"FINAL_SCHEMA_URL: {final_schema_url}\n")
        print(f"Failed to download schema from {final_schema_url}", file=sys.stderr)
        sys.exit(1)  # Script execution error, not validation error
    
    # Validate payload
    is_valid, message = validate_payload(payload, schema)
    
    # Write validation report
    status = "VALID" if is_valid else "INVALID"
    with open(args.output_report, 'w') as f:
        f.write(f"STATUS: {status}\n")
        f.write(f"MESSAGE: {message}\n")
        f.write(f"PAYLOAD_FILE: {args.payload}\n")
        f.write(f"ANALYSIS_TYPE: {analysis_type}\n")
        f.write(f"BASE_SCHEMA_URL: {args.schema_url}\n")
        f.write(f"FINAL_SCHEMA_URL: {final_schema_url}\n")
    
    # Print result
    print(f"Validation {status}: {message}")
    print(f"Analysis Type: {analysis_type}")
    print(f"Schema URL: {final_schema_url}")
    
    # Always exit with success - let Nextflow handle the exit logic based on validation result
    # The validation result is captured in the report file and exit code is determined by parsing it
    sys.exit(1 if not is_valid else 0)

if __name__ == "__main__":
    main()
