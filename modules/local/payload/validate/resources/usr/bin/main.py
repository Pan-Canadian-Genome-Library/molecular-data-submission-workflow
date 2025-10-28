#!/usr/bin/env python3

import json
import sys
import argparse
import urllib.request
import urllib.error
import re
import requests
from pathlib import Path

def download_schema(schema_url):
    """Download JSON schema from URL"""
    try:
        with urllib.request.urlopen(schema_url) as response:
            schema_content = response.read().decode('utf-8')
            return json.loads(schema_content)
    except urllib.error.URLError as e:
        print(f"Error downloading schema from {schema_url}: {e}.", file=sys.stderr)
        return None
    except json.JSONDecodeError as e:
        print(f"Error parsing schema JSON: {e}.", file=sys.stderr)
        return None

def load_payload(payload_file):
    """Load JSON payload from file"""
    try:
        with open(payload_file, 'r', encoding='utf-8') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error parsing payload JSON: {e}.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading payload file: {e}.", file=sys.stderr)
        return None
def query_clinical_validator(url,value):
        # headers={
        #         "Authorization" : "Bearer %s" % token
        # }
        #https://submission.pcgl-dev.cumulus.genomeinformatics.org/validator/category/1/entity/participant/exists?organization=EXAMPLE-CA&value=DONOR_01
        try:
                #response=requests.get(url,headers=headers)
                response=requests.get(url)
        except:
                return(False,'ERROR REACHING %s.' % (url))

        if response.status_code==404:
            return(False,"Record %s not found." % (value))
        elif response.status_code!=200:
            return(False,'ERROR w/ %s : Code %s.' % (url,response.status_code))

        if response.json()['message']=='Record found.':
            return(True,None)
        else:
            return(False,"Record %s not found." % (value))

def external_id_validate(payload,actual_schema,clinical_url):
    if "externalValidations" in actual_schema.keys():
        for validation in actual_schema['externalValidations']:
            check_url=validation['url']
            check_item=validation['jsonPath']
            study_value=payload['studyId']
            check_item_value=payload[check_item]
            ###Expected : https://submission.ingress.dev.k8s.pcgl.dev-sd4h.ca/validator/entity/experiment/field/submitter_experiment_id/exists?study={study}&value={value}
            query_url= "%s/%s" % (clinical_url,re.search(r'validator.*$',check_url)[0].replace("{study}",study_value).replace("{value}",check_item_value))
            return(query_clinical_validator(query_url,check_item_value))
    else:
        ###Nothing to check
        return True,None

    #print(payload)
    #print(actual_schema)

def validate_payload(payload, schema):
    """Validate payload against schema using jsonschema library"""
    try:
        import jsonschema
    except ImportError:
        error_msg = "jsonschema library not available in container"
        print(error_msg, file=sys.stderr)
        return False, error_msg
    
    # Extract the actual schema from the wrapper if it exists
    # The downloaded schema may be wrapped in metadata, we need the actual JSON schema
    actual_schema = schema
    if isinstance(schema, dict) and 'schema' in schema:
        actual_schema = schema['schema']
        print(f"Extracted actual schema from wrapper.", file=sys.stderr)
    else:
        print(f"Using schema as-is (no wrapper detected).", file=sys.stderr)
    
    try:
        jsonschema.validate(instance=payload, schema=actual_schema)
        return True, "Payload Validation successful."
    except jsonschema.ValidationError as e:
        return False, f"Payload Validation failed: {e.message}."
    except jsonschema.SchemaError as e:
        return False, f"Schema error: {e.message}."
    except Exception as e:
        return False, f"Validation error: {str(e)}."

def main():
    parser = argparse.ArgumentParser(description='Validate JSON payload against schema from URL')
    parser.add_argument('--payload', required=True, help='Path to JSON payload file')
    parser.add_argument('--schema-url', required=True, help='URL to download JSON schema')
    parser.add_argument('--clinical-url', required=True, help='URL to verify clinical data')
    
    args = parser.parse_args()
    
    # Load payload first to get analysisType
    payload = load_payload(args.payload)
    if payload is None:
        print("Failed to load payload file.", file=sys.stderr)
        sys.exit(1)  # Script execution error, not validation error
    
    # Extract analysisType from payload to construct final schema URL
    analysis_type = payload.get('analysisType', {}).get('name', '')
    if not analysis_type:
        print("analysisType.name not found in payload.", file=sys.stderr)
        sys.exit(1)  # Script execution error, not validation error
    
    # Construct final schema URL
    final_schema_url = f"{args.schema_url}/{analysis_type}"
    
    # Download schema
    schema = download_schema(final_schema_url)
    if schema is None:
        print(f"Failed to download schema from {final_schema_url}.", file=sys.stderr)
        sys.exit(1)  # Script execution error, not validation error
    
    # Validate payload
    messages=[]
    is_valid, message_a = validate_payload(payload, schema)
    is_valid, message_b = external_id_validate(payload, schema,args.clinical_url)

    None if message_a==None else messages.append(message_a)
    None if message_b==None else messages.append(message_b)
    
    # Print error message to stderr if validation failed
    if not is_valid:
        print("".join([msg for msg in messages]), file=sys.stderr)
    
    # Print result to stdout
    print(f"Validation {'VALID' if is_valid else 'INVALID'}: { ';'.join([msg for msg in messages if msg!=None])}")
    print(f"Analysis Type: {analysis_type}")
    print(f"Schema URL: {final_schema_url}")
    
    # Exit with validation result
    sys.exit(0 if is_valid else 1)

if __name__ == "__main__":
    main()
