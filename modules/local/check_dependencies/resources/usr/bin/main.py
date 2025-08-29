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

 Author:Edmund Su <edmund.su@oicr.on.ca>
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

def check_clinical_health(clinical_url,token):
    print("Checking Clinical URL health")
    url="%s/health" % clinical_url
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

def check_file_manager_health(file_manager,token):
    print("Checking File Manager health")
    url="%s/isAlive" % file_manager
    # headers={
    #     "Authorization" : "Bearer %s" % token
    # }

    try:
        response=requests.get(url)
    except:
        raise ValueError('ERROR REACHING %s' % (url))

    if response.status_code!=200:
        raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
        exit(1)

def check_clinical_study(clinical_url,study_id,token):
    print("Checking Study in Clinical")
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

def check_file_manager_study(file_manager_url,study_id,token):
    print("Checking Study in File Manager")
    url="%s/studies/%s" % (file_manager_url,study_id)
    #headers={
    #    "Authorization" : "Bearer %s" % token
    #}
    try:
        response=requests.get(url)
    except:
        raise ValueError('ERROR REACHING %s' % (url))

    if response.status_code!=200:
        raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
        exit(1)

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

def check_analysis_types(file_manager_url,study_id,token):
    analysis_types=[]
    required_analysis_fields={}

    print("Retrieving analysis Types")
    # headers={
    #     "Authorization" : "Bearer %s" % token
    # }
    limit=20
    url="%s/schemas?hideSchema=true&limit=%s&offset=0&unrenderedOnly=false" % (file_manager_url,str(limit))

    try:
        response=requests.get(url)
    except:
        raise ValueError('ERROR REACHING %s' % (url))

    if response.status_code!=200:
        raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
        exit(1)

    if response.json()['count']<limit:
        analysis_types = list(set([analysis_type['name'] for analysis_type in response.json()['resultSet']]))
    else:
        total=response.json()['count']
        for offset in range(0,total,limit):
            url="%s/schemas?hideSchema=true&limit=%s&offset=%s&unrenderedOnly=false" % (file_manager_url,str(limit),str(offset))
            try:
                response=requests.get(url,headers=headers)
            except:
                raise ValueError('ERROR REACHING %s' % (url))

            if response.status_code!=200:
                raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
                exit(1)

            for analysis_type in response.json()['resultSet']:
                if analysis_type['name'] not in analysis_types:
                    analysis_types.append(analysis_type['name'])


    for analysis_type in analysis_types:
        required_analysis_fields[analysis_type]={}
        required_analysis_fields[analysis_type]['fields']={}
        required_analysis_fields[analysis_type]['dataTypes']=[]
        print("Retrieving info for %s" % analysis_type)
        url="%s/schemas/%s?unrenderedOnly=false" % (file_manager_url,analysis_type)
        try:
            response=requests.get(url)
        except:
            raise ValueError('ERROR REACHING %s' % (url))

        if response.status_code!=200:
            raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
            exit(1)

        for required in response.json()['schema']['required']:
            if response.json()['schema']['properties'][required].get('type'):
                if response.json()['schema']['properties'][required]["type"]=="string":
                    required_analysis_fields[analysis_type]['fields'][required]=[required]
                elif response.json()['schema']['properties'][required]["type"]=="array":
                    definitions=response.json()['schema']['properties'][required]['items']["$ref"].split("/")[1]
                    parent=response.json()['schema']['properties'][required]['items']["$ref"].split("/")[2]
                    root=response.json()['schema']['properties'][required]['items']["$ref"].split("/")[3]
                    required_analysis_fields[analysis_type]['fields'][required]=response.json()['schema'][definitions][parent][root]['required']
                elif response.json()['schema']['properties'][required]["type"]=="object":
                    required_analysis_fields[analysis_type]['fields'][required]=response.json()['schema']['properties'][required]['required']
            elif response.json()['schema']['properties'][required].get('allOf'):
                ###working off a very naive assumption for now on how allOf is structured relative to analysisType. Does not handle future cases
                definitions=response.json()['schema']['properties'][required]['allOf'][0]["$ref"].split("/")[1]
                parent=response.json()['schema']['properties'][required]['allOf'][0]["$ref"].split("/")[2]
                required_analysis_fields[analysis_type]['fields'][required]=response.json()['schema'][definitions][parent]['required']
            else:
                raise ValueError('ERROR w/ schema %s : missing support for types %s' % (analysis_type,",".join(response.json()['schema'][required].keys())))
                exit(1)
        
        if len(response.json()['fileTypes'])>0:
            required_analysis_fields[analysis_type]['dataTypes']=response.json()['fileTypes']
            required_analysis_fields[analysis_type]['externalValidations']=[ external for external in response.json()['externalValidations'] if external.get("url") and external.get("jsonPath")]
        
    return(required_analysis_fields)

def retrieve_clinical_schema(clinical_url,category_id,token):
    print("Retrieving Schema")
    url="%s/dictionary/category/%s" % (clinical_url,category_id)
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

    return(response.json()['schemas'])

def generate_relational_mapping(schema,check_entities):
    print("Determining relational mapping")
    
    relational_mappings={}

    for key in check_entities:
        relational_mappings[key.lower()]={}
        relational_mappings[key]['primary']=[]
        relational_mappings[key]['foreign']=[]

    for entity in schema:
    ####Check unique Id
        if entity['name'].lower() in check_entities:
            ####Check unique Id
            for field in entity['fields']:
                if field.get('unique'):
                    if field['unique']==True:
                        relational_mappings.get(entity['name']).get('primary').append(field['name'].lower())
            ####Check foreignkeys
            if entity.get("restrictions"):
                if entity.get("restrictions").get("foreignKey"):
                    for foreignKey in entity.get("restrictions").get("foreignKey"):
                        for mappings in foreignKey['mappings']:
                            relational_mappings[entity['name']]['foreign'].append(
                                {
                                    "entity":foreignKey['schema'].lower(),
                                    "foreign":mappings['local'].lower()
                                }
                            )
    return(relational_mappings)     
def check_minimum_columns(metadata_file,cols):
    print("Checking minimum columns for %s" % metadata_file)
    data=pd.read_csv(metadata_file,sep='\t')

    for col in cols :
        if col not in data.columns.values.tolist():
            raise ValueError('col %s is missing in file %s' % (col,metadata_file))
            exit(1)

def update_relational_mapping(relational_mapping,analysis_types):
    print("Adding additional mapping for analysis,workflow and files")
    relational_mapping['files']={
        "primary":["fileName"],
        "foreign":[{
            "foreign" : "submitter_analysis_id",
            "entity" : "analysis"
        }]
    }

    relational_mapping['analysis']={"analysisTypes":{}}

    for schema in analysis_types:
        #print(schema)
        if analysis_types[schema].get('externalValidations'):
            if len(analysis_types[schema]['externalValidations'])==0:
                continue
            if analysis_types[schema]['externalValidations'][0]['url'].endswith("health"):
                analysis_types[schema]['externalValidations'][0]['url']="https://submission.pcgl-dev.cumulus.genomeinformatics.org/validator/category/1/entity/experiment/exists?organization={studyId}&value={value}"
            
            relational_mapping['analysis']['analysisTypes'][schema]={
                "primary":["submitter_analysis_id"],
                "foreign":{
                    "foreign":analysis_types[schema]['externalValidations'][0]['jsonPath'],
                    "entity":re.findall(r'(?<=entity\/)([^\/]+)(?=\/exists)',analysis_types[schema]['externalValidations'][0]['url'])[0]
                }
            }
    
    relational_mapping['workflow']={"analysisTypes":{}}

    for schema in analysis_types:
        if analysis_types[schema].get("fields").get("workflow"):
            relational_mapping['workflow']['analysisTypes'][schema]={}

    

def main(args):
    if args.file_metadata: print("input:",args.file_metadata)
    if args.analysis_metadata: print("input:",args.analysis_metadata)
    if args.workflow_metadata: print("input:",args.workflow_metadata)
    if args.sample_metadata: print("input:",args.sample_metadata)
    if args.specimen_metadata: print("input:",args.specimen_metadata)
    if args.experiment_metadata: print("input:",args.experiment_metadata)
    if args.read_group_metadata: print("input:",args.read_group_metadata)
    if args.clinical_url: print("input:",args.clinical_url)
    if args.file_manager_url: print("input:",args.file_manager_url)
    if args.study_id: print("input:",args.study_id)
    if args.token: print("input:",args.token)

    ###Preflight checks
    check_clinical_health(
        args.clinical_url,
        args.token
        )

    check_file_manager_health(
        args.file_manager_url,
        args.token
        )

    ###R1b - The pipeline shall query the study service to check for registered study. Study registration is assumed pre-existing; failure here should halt the pipeline.
    check_clinical_study(
           args.clinical_url,
           args.study_id,
           args.token
           )

    check_file_manager_study(
           args.file_manager_url,
           args.study_id,
           args.token
           )

    ###Retrieve requirements
    category_id=retrieve_category_id(
        args.clinical_url,
        args.study_id,
        args.token
        )
        
    clinical_schema=retrieve_clinical_schema(
        args.clinical_url,
        category_id,
        args.token
        )

    analysis_types=check_analysis_types(
           args.file_manager_url,
           args.study_id,
           args.token
           )

    ###Retrieve foreign dependencies
    relational_mapping=generate_relational_mapping(
        clinical_schema,
        ["experiment","read_group","sample","specimen"]
    )

    ###Check minimum columns to infer relationships
    for metadata,key in zip(
        [args.sample_metadata,args.specimen_metadata,args.experiment_metadata,args.read_group_metadata],
        ["sample","specimen","experiment","read_group"],
    ):
        if metadata:
            check_minimum_columns(
                metadata,
                relational_mapping.get(key).get('primary')+ [foreign_key.get('foreign') for foreign_key in relational_mapping.get(key).get('foreign')]
            )

    if args.workflow_metadata:
        check_minimum_columns(
            args.workflow_metadata,
             ["submitter_workflow_id"]
        )
    check_minimum_columns(
        args.analysis_metadata,
         ["analysisType","submitter_analysis_id"]
    )

    check_minimum_columns(
        args.file_metadata,
        ["fileName","dataType",'fileMd5sum']
    )

    update_relational_mapping(relational_mapping,analysis_types)

    with open('relational_mapping.json', 'w') as f:
        json.dump(relational_mapping, f)

    with open("analysis_types.json","w") as f:
        json.dump(analysis_types, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool: Check Clinical dependencies')
    parser.add_argument("-fi", "--file_metadata", dest="file_metadata", required=True, help="file metadata tsv")
    parser.add_argument("-an", "--analysis_metadata", dest="analysis_metadata", required=True, help="analysis metadata tsv")
    parser.add_argument("-wo", "--workflow_metadata", default=False, dest="workflow_metadata", required=False, help="workflow metadata tsv")
    parser.add_argument("-sa", "--sample_metadata", default=False, dest="sample_metadata", required=False, help="sample metadata tsv")
    parser.add_argument("-sp", "--specimen_metadata", default=False, dest="specimen_metadata", required=False, help="specimen metadata tsv")
    parser.add_argument("-ex", "--experiment_metadata", default=False, dest="experiment_metadata", required=False, help="experiment metadata tsv")
    parser.add_argument("-rg", "--read_group_metadata", default=False, dest="read_group_metadata", required=False, help="read_group metadata tsv")
    parser.add_argument("-cu", "--clinical_url", dest="clinical_url", required=True, help="Clinical URL")
    parser.add_argument("-fm", "--file_manager_url", dest="file_manager_url", required=True, help="File Manager URL")
    parser.add_argument("-si", "--study_id", dest="study_id", required=True, help="study_id")
    parser.add_argument("-t", "--token", dest="token", required=True, help="token")

    args = parser.parse_args()

    main(args)