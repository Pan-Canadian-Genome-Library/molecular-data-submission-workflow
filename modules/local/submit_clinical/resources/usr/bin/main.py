#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import requests
import numpy as np
import re
import jsonschema
import json
import argparse
import os
import glob
from typing import Dict, List
import time

def retrieve_category_id(clinical_url,study_id,token):
    print("Retrieve Category ID")
    url="%s/study/%s" % (clinical_url,study_id)

    headers={
            "Authorization" : "Bearer %s" % token
    }
    try:
            response=requests.get(url,headers=headers)
    except:
            raise ValueError('ERROR REACHING %s' % (url))

    if response.status_code!=200:
            raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
            exit(1)

    if response.json().get('categoryId'):
        if response.json().get('categoryId')!=None:
            return(str(response.json().get('categoryId')))
        else:
            raise ValueError('ERROR w/ %s : %s study\'s corresponding schema was not found ' % (url,study_id))
    else:
        raise ValueError('ERROR w/ %s : %s study\'s corresponding schema was not found ' % (url,study_id))

    for cat_id in categories:
        if study_id.lower() in cat_id['name'] or study_id.upper() in cat_id['name']:
            return(str(cat_id["id"]))

    for cat_id in categories:
        if "prod_pcgl_schema" in cat_id['name']:
            return(str(cat_id["id"]))

    raise ValueError('ERROR w/ %s : %s study\'s corresponding schema was not found ' % (url,study_id))


def check_existing_submission(category_id,clinical_url,study_id,token):
    print("Verifying existance of existing submission")
    #https://submission.pcgl-dev.cumulus.genomeinformatics.org/submission/category/1?onlyActive=true&organization=EXAMPLE-CA
    url="%s/submission/category/%s?onlyActive=true&organization=%s" % (clinical_url,str(category_id),study_id)
    headers={
            "Authorization" : "Bearer %s" % token
    } 
    try:
            response=requests.get(url,headers=headers)
    except:
        comments=[]
        comments.append('ERROR REACHING %s' % (url))

        comments.append(response.json().get('error')) if response.json().get('error') else comments
        comments.append("message : %s " % response.json().get('message')) if response.json().get('message') else comments
        
        raise ValueError("\n".join(comments))
        exit(1)

    
    if response.status_code==404:
        print("No existing submission found")
        return(True)

    if response.status_code!=200:
        comments=[]
        comments.append('ERROR w/ %s : Code %s' % (url,response.status_code))

        comments.append(response.json().get('error')) if response.json().get('error') else comments
        comments.append("message : %s " % response.json().get('message')) if response.json().get('message') else comments
        
        raise ValueError("\n".join(comments))
        exit(1)

    if response.json()['pagination']['totalRecords']>0:
        for submisison_id in [record['id'] for record in response.json()['records']]:
                delete_existing_submission(clinical_url,submisison_id,token)
        check_existing_submission(category_id,clinical_url,study_id,token)
    
    return(True)

def delete_existing_submission(clinical_url,submisison_id,token):
    print("Deleting previous submission to enable new submission")
    #https://submission.pcgl-dev.cumulus.genomeinformatics.org/submission/30
    url="%s/submission/%s" % (clinical_url,submisison_id)
    headers={
            "Authorization" : "Bearer %s" % token
    }

    try:
            response=requests.delete(url,headers=headers)
    except:
        comments=[]
        comments.append('ERROR REACHING %s' % (url))

        comments.append(response.json().get('error')) if response.json().get('error') else comments
        comments.append("message : %s " % response.json().get('message')) if response.json().get('message') else comments
        
        raise ValueError("\n".join(comments))
        exit(1)

    if response.status_code!=200:
        comments=[]
        comments.append('ERROR w/ %s : Code %s' % (url,response.status_code))

        comments.append(response.json().get('error')) if response.json().get('error') else comments
        comments.append("message : %s " % response.json().get('message')) if response.json().get('message') else comments
        
        raise ValueError("\n".join(comments))
        exit(1)

    return(True)

def read_data(sample_metadata,specimen_metadata,experiment_metadata,read_group_metadata):
    print("Ingesting biospecimen entities")
    analysis={}
    analysis['comments']=[]
    for metadata,entity in zip([sample_metadata,specimen_metadata,experiment_metadata,read_group_metadata],["sample","specimen","experiment","read_group"]):
        if metadata:
            tmp=pd.read_csv(metadata,sep='\t',index_col=False)
            analysis[entity]=tmp.copy()

            for col in analysis[entity].columns.values.tolist():
                if len(analysis[entity])>1:
                    if False not in pd.isna(analysis[entity][col].values.tolist()):
                        analysis[entity].drop(col,axis=1,inplace=True)
                else:
                    if pd.isna(analysis[entity].loc[0,col]):
                        analysis[entity].drop(col,axis=1,inplace=True)
    return(analysis)

def query_clinical_validator(url,token):
        headers={
                "Authorization" : "Bearer %s" % token
        }


        #https://submission.pcgl-dev.cumulus.genomeinformatics.org/validator/category/1/entity/participant/exists?organization=EXAMPLE-CA&value=DONOR_01
        try:
                response=requests.get(url,headers=headers)
        except:
                raise ValueError('ERROR REACHING %s' % (url))

        if response.status_code!=200 and response.status_code!=404:
                raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
                exit(1)

        if response.json()['message']=='Record found.':
            return(True)
        else:
            return(False)


def query_registered_data(token,clinical_url,category_id,entity,study_id,primary_key,val):
    headers={
            "Authorization" : "Bearer %s" % token
    }
    url="%s/data/category/%s/organization/%s/query?entityName=%s" % (clinical_url,category_id,study_id,entity)
    payload={
            "op": "and",
            "content": [
                    {
                            "op": "in",
                            "content": {
                                    "fieldName": primary_key,
                                    "value": [val]
                                    }
                            }
                    ]
            }
    print(url)
    print(payload)
    try:
            response=requests.post(url,json=payload,headers=headers)
    except:
            raise ValueError('ERROR REACHING %s' % (url))

    if response.status_code!=200:
            raise ValueError('ERROR w/ %s : Code %s - %s' % (url,response.status_code,response.text))
            exit(1)
    return(response)  

def verify_registered_data(token,clinical_url,category_id,entity,study_id,primary_key,ind,analysis):
    print("Verifying if submitted %s data record %s is consistent" % (entity,analysis[entity].loc[ind,primary_key]))
    status=[]
    comments=[]

    response=query_registered_data(token,clinical_url,category_id,entity,study_id,primary_key,analysis[entity].loc[ind,primary_key])

    for col in analysis[entity].columns.values.tolist():
        if response.json()['records'][0]['data'].get(col):
            valA=response.json()['records'][0]['data'][col]

            valB=analysis[entity].loc[ind,col]

            if valA!=valB:
                comments.append("Field '%s' is not consistent for record %s in entity %s. Specified - %s vs Comitted - %s" % (col,analysis[entity].loc[ind,primary_key],entity,valA,valB))
                status.append(False)
        else:
            comments.append("New Field '%s' detected for existing record %s in entity %s. Contact PCGL Admin to perform this update seperately" % (col,analysis[entity].loc[ind,primary_key],entity))
            status.append(False)

    return(comments,status)
    

def check_registered_entities(analysis,clinical_url,category_id,study_id,relational_mapping,token):
    print("Checking registration status for biospecimen entities")
    usable={} ### Tracks entity usage. If entity record has no changes and is redundant, we avoid resubmitting
    for entity in analysis.keys():
        if entity=='comments':
            continue
        for primary_key in relational_mapping[entity]['primary']:
            ###Should be 1:1 aside from read_group, otherwise if one read_group fails all will be flagged
            for ind in analysis[entity].index.values.tolist():
                url="%s/validator/entity/%s/field/%s/exists?study=%s&value=%s" % (clinical_url,entity,primary_key,study_id,analysis[entity].loc[ind,primary_key])
                print(url)
                #https://submission.pcgl-dev.cumulus.genomeinformatics.org/validator/entity/experiment/field/submitter_experiment_id/exists?study={study}&value={value}
                if query_clinical_validator(url,token):
                    comments,redundant=verify_registered_data(token,clinical_url,category_id,entity,study_id,primary_key,ind,analysis)

                    ###Regardless of record matching to committed data, we will not submit
                    usable[entity]=False

                    for comment in comments:
                        analysis['comments'].append(comment)
                else:
                    usable[entity]=True

    return usable


def rename_input(output_directory,sample_metadata,specimen_metadata,experiment_metadata,read_group_metadata,usability):
    print("Renaming files to fit lyric")
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    if not os.path.isdir("%s/submit" % (output_directory)):
        os.makedirs("%s/submit" % (output_directory))

    for metadata,entity in zip([sample_metadata,specimen_metadata,experiment_metadata,read_group_metadata],["sample","specimen","experiment","read_group"]):
        ###We only submit entities that present and not redundant (i.e. submitted before)
        if metadata and usability[entity]:
            tmp=pd.read_csv(metadata,sep='\t',index_col=False)
            rename=os.path.basename(entity).lower().capitalize()
            print("Renaming %s to %s" % (metadata,rename))
            tmp.to_csv("%s/submit/%s.tsv" % (output_directory,rename),sep='\t',index=False)
    return(True)

def submit_clinical(clinical_url,category_id,study_id,output_directory,token):
    print("Using Clinical-Submission endpoint to submit data")
    url='%s/submission/category/%s/data' % (clinical_url,str(category_id))
    headers={
            "Authorization" : "Bearer %s" % token,
            'accept': 'application/json'
    } 

    files=[]
    for file in glob.iglob("%s/submit/*.tsv" % (output_directory)):
        files.append(
            (
                'files',
                (
                    file.replace(output_directory+"/",""),
                    open(file,'rb'),
                    'text/tab-separated-values'
                )
            )
        )


    try:
            response = requests.post(url, headers=headers, files=files,data={"organization":study_id})
            print(files)
            print(headers)
    except:
        comments=[]
        comments.append('ERROR REACHING %s' % (url))

        comments.append("error : %s" % response.json().get('error')) if response.json().get('error') else comments
        comments.append("message : %s " % response.json().get('message')) if response.json().get('message') else comments
        
        raise ValueError(";\n".join(comments))
        exit(1)

    if response.status_code!=200:
        comments=[]
        comments.append('ERROR w/ %s : Code %s' % (url,response.status_code))

        comments.append("error : %s" % response.json().get('error')) if response.json().get('error') else comments
        comments.append("message : %s " % response.json().get('message')) if response.json().get('message') else comments
        print(response.json())
        raise ValueError(";\n".join(comments))
        exit(1)

    return(str(response.json()['submissionId']))

def check_submission_status(
    clinical_url: str,
    submission_id: str,
    token: str,
    stage: str = 'validation',
    poll_interval: int = 5,
    max_wait: int = 300
) -> bool:
    """
    Poll submission status until it reaches the expected terminal state.

    Validation stage: polls until VALID (waits on OPEN).
    Commit stage:     polls until COMMITTED (waits on VALID while server commits).
    Raises on INVALID in either stage.

    Status enum: OPEN, VALIDATING, VALID, INVALID, CLOSED, COMMITTING, COMMITTED

    Args:
        clinical_url: Base URL for clinical API
        submission_id: Submission ID to check
        token: Authentication token
        stage: Either 'validation' or 'commit'
        poll_interval: Seconds between status polls
        max_wait: Maximum total seconds to wait before raising

    Returns:
        True when the expected terminal status is reached

    Raises:
        ValueError: If status is INVALID or CLOSED, an unexpected state, or max_wait exceeded
    """
    stage_msg = "Validating" if stage == 'validation' else "Verifying committed"
    print(f"{stage_msg} submission: {submission_id}")

    elapsed = 0
    while elapsed < max_wait:
        response = api_request('GET', f"{clinical_url}/submission/{submission_id}", token)
        status = response.json()['status']

        if status in ['INVALID', 'CLOSED']:
            response_details = api_request('GET', f"{clinical_url}/submission/{submission_id}/details", token)
            errors = parse_errors_detail(response_details.json())
            stage_error = "Validation" if stage == 'validation' else "Commit"
            raise ValueError(f"{stage_error} failed with errors:\n" + "\n".join(errors))

        if stage == 'validation':
            if status == 'VALID':
                print("Validation successful")
                return True
            if status in ['OPEN', 'VALIDATING']:
                print(f"Status: {status} — waiting... ({elapsed}s elapsed)")
            else:
                raise ValueError(
                    f"Submission {submission_id} reached unexpected status '{status}' "
                    f"during validation stage."
                )

        elif stage == 'commit':
            if status == 'COMMITTED':
                print("Data successfully committed to database")
                return True
            if status in ['VALID', 'COMMITTING']:
                # Server is processing the commit
                print(f"Status: {status} — commit in progress... ({elapsed}s elapsed)")
            else:
                raise ValueError(
                    f"Submission {submission_id} reached unexpected status '{status}' "
                    f"during commit stage."
                )

        time.sleep(poll_interval)
        elapsed += poll_interval

    raise ValueError(
        f"Timed out after {max_wait}s waiting for submission {submission_id} "
        f"to reach {'VALID' if stage == 'validation' else 'COMMITTED'} status."
    )

def parse_errors_detail(response_data: dict) -> List[str]:
    """
    Parse errors details from API response.
    
    Args:
        response_data: JSON response data
        
    Returns:
        List of formatted error messages
    """
    errors = []
    insert_errors = response_data.get('errors', [])
    
    for error in insert_errors:
        field_name = error.get('fieldName', 'N/A')
        reason = error.get('reason', 'Unknown')
        field_value = error.get('fieldValue', '')
        
        error_parts = [f"Reason: {reason}", f"Field: {field_name}"]
        if field_value:
            error_parts.append(f"Value: {field_value}")
        errors.append("  " + ", ".join(error_parts))
    
    return errors

def api_request(
    method: str,
    url: str,
    token: str,
    timeout: int = 30,
    **kwargs
) -> requests.Response:
    """
    Make authenticated API request with error handling.
    
    Args:
        method: HTTP method (GET, POST, DELETE)
        url: API endpoint URL
        token: Authentication token
        timeout: Request timeout in seconds
        **kwargs: Additional arguments for requests
        
    Returns:
        Response object
        
    Raises:
        ValueError: If request fails
    """
    headers = kwargs.pop('headers', {})
    headers.setdefault('Authorization', f'Bearer {token}')
    headers.setdefault('accept', 'application/json')
    
    try:
        response = requests.request(method, url, headers=headers, timeout=timeout, **kwargs)
    except requests.exceptions.RequestException as e:
        raise ValueError(f'ERROR reaching {url}: {e}')
    
    if response.status_code not in [200, 404]:
        # handle_api_error
        comments = [f'ERROR with {url}: Code {response.status_code}']
    
        try:
            error_data = response.json()
            if error_data.get('error'):
                comments.append(f"Error: {error_data['error']}")
            if error_data.get('message'):
                comments.append(f"Message: {error_data['message']}")
        except Exception:
            comments.append(f"Response: {response.text}")
        
        raise ValueError("\n".join(comments))
    
    return response

def commit_clinical(clinical_url,category_id,submission_id,token):
    print("Written submission contents to database")
    url="%s/submission/category/%s/commit/%s" % (clinical_url,category_id,submission_id)
    
    headers={
            "Authorization" : "Bearer %s" % token,
            'accept': 'application/json'
    } 
    try:
        response = requests.post(url, headers=headers)
    except:
        comments=[]
        comments.append('ERROR REACHING %s' % (url))

        comments.append(response.json().get('error')) if response.json().get('error') else comments
        comments.append("message : %s " % response.json().get('message')) if response.json().get('message') else comments
        
        raise ValueError("\n".join(comments))
        exit(1)

    if response.status_code!=200:
        comments=[]
        comments.append('ERROR w/ %s : Code %s' % (url,response.status_code))

        comments.append(response.json().get('error')) if response.json().get('error') else comments
        comments.append("message : %s " % response.json().get('message')) if response.json().get('message') else comments
        
        raise ValueError("\n".join(comments))
        exit(1)

    return(True)


def return_submitted_data(
    category_id,
    clinical_url,
    token,
    analysis,
    output_directory
):
    print("Check if submitted clinical data is needed for downstream molecular validation")

    data=pd.read_csv(analysis,sep='\t',index_col=False)
    output={}

    for ind in data.index.values.tolist():
        if data.loc[ind,"analysisType"]=='sequenceExperiment':
            print("sequenceExperiment confirmed. Clinical data is needed for downstream molecular validation")
            ###Currently does not account for pagination
            for primary_key,entity in zip(["submitter_experiment_id"],['read_group']):
                response=query_registered_data(
                    token,
                    clinical_url,
                    category_id,
                    entity,
                    data.loc[ind,"studyId"],
                    primary_key,
                    data.loc[ind,primary_key]
                    )
                
                print(response.text)
                if response.status_code!=404:
                    output[entity]=pd.DataFrame()

                for record in response.json()['records']:
                    ind=len(output[entity])
                    for key in record['data'].keys():
                        output[entity].loc[ind,key]=record['data'][key]
        else:
            print("Verifying Dependency for analysis record %s is met" % data.loc[ind,"submitter_analysis_id"])
            ###Currently does not account for pagination
            for primary_key,entity in zip(["submitter_experiment_id"],['experiment']):
                response=query_registered_data(
                    token,
                    clinical_url,
                    category_id,
                    entity,
                    data.loc[ind,"studyId"],
                    primary_key,
                    data.loc[ind,primary_key]
                    )
                  

    if len(output.keys())>0:
    
        if not os.path.isdir(output_directory):
            os.makedirs(output_directory)

        if not os.path.isdir("%s/retrieved" % (output_directory)):
            os.makedirs("%s/retrieved" % (output_directory))

        for entity in output.keys():
            print("Outputing %s as %s" % (entity,"%s/retrieved/%s.tsv" % (output_directory,entity)))
            output[entity].to_csv("%s/retrieved/%s.tsv" % (output_directory,entity),sep='\t',index=False)
    else:
        print("Nothing to retrieve")

def main(args):
    if args.analysis_metadata: print("analysis_metadata:",args.analysis_metadata)
    if args.sample_metadata: print("sample_metadata:",args.sample_metadata)
    if args.specimen_metadata: print("specimen_metadata:",args.specimen_metadata)
    if args.experiment_metadata: print("experiment_metadata:",args.experiment_metadata)
    if args.read_group_metadata: print("read_group_metadata:",args.read_group_metadata)
    if args.clinical_url: print("clinical_url:",args.clinical_url)
    if args.study_id: print("study_id:",args.study_id)
    if args.output_directory: print("output_directory:",args.output_directory)
    if args.relational_mapping: print("relational_mapping:",args.relational_mapping)

    category_id=retrieve_category_id(
        args.clinical_url,
        args.study_id,
        args.token
        )
    print(category_id)

    check_existing_submission(
        category_id,
        args.clinical_url,
        args.study_id,
        args.token
    )

    with open(args.relational_mapping, 'r') as file:
        relational_mapping = json.load(file)

    analysis=read_data(
        args.sample_metadata,
        args.specimen_metadata,
        args.experiment_metadata,
        args.read_group_metadata
    )

    usable=check_registered_entities(analysis,args.clinical_url,category_id,args.study_id,relational_mapping,args.token) ###Returns an array of booleans, bool status reflects if any delta in columns or entity record was previously submitted with no changes (redundant submission)

    rename_input(
        args.output_directory,
        args.sample_metadata,
        args.specimen_metadata,
        args.experiment_metadata,
        args.read_group_metadata,
        usable
    )
    
    if len(analysis['comments'])>0:
        raise ValueError(".\n".join(analysis['comments']))

    if len([file for file in glob.iglob(args.output_directory+"/submit/*.tsv")])>0:
        submission_id=submit_clinical(
            args.clinical_url,
            category_id,
            args.study_id,
            args.output_directory,
            args.token
        )

        print(submission_id)

        check_submission_status(
            args.clinical_url,
            submission_id,
            args.token,
            stage='validation'
        )

        commit_clinical(
            args.clinical_url,
            category_id,
            submission_id,
            args.token
        )

        check_submission_status(
            args.clinical_url,
            submission_id,
            args.token,
            stage='commit'
        )

    else:
        print("No data to submit")

    return_submitted_data(
        category_id,
        args.clinical_url,
        args.token,
        args.analysis_metadata,
        args.output_directory
    )



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool: Check Clinical dependencies')
    parser.add_argument("-an", "--analysis_metadata", dest="analysis_metadata", required=True, help="analysis metadata tsv")
    parser.add_argument("-sa", "--sample_metadata", default=False, dest="sample_metadata", required=False, help="sample metadata tsv")
    parser.add_argument("-sp", "--specimen_metadata", default=False, dest="specimen_metadata", required=False, help="specimen metadata tsv")
    parser.add_argument("-ex", "--experiment_metadata", default=False, dest="experiment_metadata", required=False, help="experiment metadata tsv")
    parser.add_argument("-rg", "--read_group_metadata", default=False, dest="read_group_metadata", required=False, help="read_group metadata tsv")
    parser.add_argument("-cu", "--clinical_url", dest="clinical_url", required=True, help="Clinical URL")
    parser.add_argument("-si", "--study_id", dest="study_id", required=True, help="study_id")
    parser.add_argument("-t", "--token", dest="token", required=True, help="token")
    parser.add_argument("-od", "--output-directory", dest="output_directory", required=False, help="output directory where entity files are saved by analysis",default="output")
    parser.add_argument("-rm", "--relational_mapping", dest="relational_mapping", required=True, help="relational mapping json")
    
    args = parser.parse_args()

    main(args)