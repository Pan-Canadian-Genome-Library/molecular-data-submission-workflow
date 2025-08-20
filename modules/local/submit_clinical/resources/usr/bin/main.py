#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Copyright (c) 2019, Ontario Institute for Cancer Research (OICR).

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published
 by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

 Author:Guanqiao Feng <gfeng@oicr.on.ca>
        Junjun Zhang <junjun.zhang@oicr.on.ca>
        Linda Xiang <linda.xiang@oicr.on.ca>
 """

import pandas as pd
import requests
import numpy as np
import re
import jsonschema
import json
import argparse
import os
import glob

def retrieve_category_id(clinical_url,study_id,token):
        print("Retrieve Category ID")
        url="%s/category" % (clinical_url)
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

        categories=response.json()

        for cat_id in categories:
                if study_id.lower() in cat_id['name']:
                        return(str(cat_id["id"]))

        for cat_id in categories:
                if "prod_pcgl_schema" in cat_id['name']:
                        return(str(cat_id["id"]))

        raise ValueError('ERROR w/ %s : %s study\'s corresponding schema was not found ' % (url,study_id))


def check_existing_submission(category_id,clinical_url,study_id,token):
    print("check_existing_submission")
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
        print("GOOD NEWS NOTHING FOUND")
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
    print("delete_existing_submission")
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

def rename_input(output_directory,sample_metadata,specimen_metadata,experiment_metadata,read_group_metadata):
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    for entity in [sample_metadata,specimen_metadata,experiment_metadata,read_group_metadata]:
        if entity:
            tmp=pd.read_csv(entity,sep='\t')
            rename=os.path.basename(entity).lower().capitalize()
            print(rename)
            tmp.to_csv("%s/%s.tsv" % (output_directory,rename),sep='\t',index=False)
    return(True)

def submit_clinical(clinical_url,category_id,study_id,output_directory,token):
    print("submit_clinical")
    url='%s/submission/category/%s/data' % (clinical_url,str(category_id))
    headers={
            "Authorization" : "Bearer %s" % token,
            'accept': 'application/json'
    } 

    files=[]
    for file in glob.iglob(output_directory+"/*.tsv"):
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
    print(response.json())
    return(str(response.json()['submissionId']))

def validation_status(clinical_url,submission_id,token):
    print("validation_status")
    url="%s/submission/%s" % (clinical_url,submission_id)
    headers={
            "Authorization" : "Bearer %s" % token,
            'accept': 'application/json'
    } 

    try:
            response = requests.get(url, headers=headers)
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

    if response.json()['status']=='INVALID':
        comments=[]

        for entity in response.json()['errors']['inserts'].keys():
            for error in response.json()['errors']['inserts'][entity]:
                if error.get('fieldValue'):
                    comments.append("Entity - %s:%s, field:%s, value:%s" % (entity,error['reason'],error['fieldName'],error['fieldValue']))
                else:
                    comments.append("Entity - %s:%s, field:%s" % (entity,error['reason'],error['fieldName']))
        raise ValueError("\n".join(comments))
        exit(1)
    return(True)

def commit_clinical(clinical_url,category_id,submission_id,token):
    print("commit_clinical")
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

def committed_status(clinical_url,submission_id,token):
    print("committed_status")
    url="%s/submission/%s" % (clinical_url,submission_id)
    
    headers={
            "Authorization" : "Bearer %s" % token,
            'accept': 'application/json'
    } 

    try:
            response = requests.get(url, headers=headers)
    except:
            raise ValueError('ERROR REACHING %s' % (url))

    if response.status_code!=200:
            raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
            exit(1)

    if response.json()['status']=='INVALID':
        comments=[]

        for entity in response.json()['errors']['inserts'].keys():
            for error in response.json()['errors']['inserts'][entity]:
                comments.append("%s %s %s %s" % (entity,error['reason'],error['fieldName'],error['fieldValue']))

        raise ValueError("\n".join(comments))
        exit(1)
    return(True)

def main(args):
    if args.sample_metadata: print("input:",args.sample_metadata)
    if args.specimen_metadata: print("input:",args.specimen_metadata)
    if args.experiment_metadata: print("input:",args.experiment_metadata)
    if args.read_group_metadata: print("input:",args.read_group_metadata)
    if args.clinical_url: print("input:",args.clinical_url)
    if args.study_id: print("input:",args.study_id)
    if args.token: print("input:",args.token)
    if args.output_directory: print("input:",args.output_directory)

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

    rename_input(
        args.output_directory,
        args.sample_metadata,
        args.specimen_metadata,
        args.experiment_metadata,
        args.read_group_metadata
    )

    submission_id=submit_clinical(
        args.clinical_url,
        category_id,
        args.study_id,
        args.output_directory,
        args.token
    )

    print(submission_id)

    validation_status(
        args.clinical_url,
        submission_id,
        args.token
    )

    commit_clinical(
        args.clinical_url,
        category_id,
        submission_id,
        args.token
    )

    committed_status(        
        args.clinical_url,
        submission_id,
        args.token
    )






if __name__ == "__main__":
        parser = argparse.ArgumentParser(description='Tool: Check Clinical dependencies')
        parser.add_argument("-sa", "--sample_metadata", default=False, dest="sample_metadata", required=False, help="sample metadata tsv")
        parser.add_argument("-sp", "--specimen_metadata", default=False, dest="specimen_metadata", required=False, help="specimen metadata tsv")
        parser.add_argument("-ex", "--experiment_metadata", default=False, dest="experiment_metadata", required=False, help="experiment metadata tsv")
        parser.add_argument("-rg", "--read_group_metadata", default=False, dest="read_group_metadata", required=False, help="read_group metadata tsv")
        parser.add_argument("-cu", "--clinical_url", dest="clinical_url", required=True, help="Clinical URL")
        parser.add_argument("-si", "--study_id", dest="study_id", required=True, help="study_id")
        parser.add_argument("-t", "--token", dest="token", required=True, help="token")
        parser.add_argument("-od", "--output-directory", dest="output_directory", required=False, help="output directory where entity files are saved by analysis",default="output")
        
        args = parser.parse_args()

        main(args)