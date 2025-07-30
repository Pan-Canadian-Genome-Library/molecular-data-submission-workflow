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
                url="%s/category/%s" % (clinical_url,str(cat_id['id']))

                try:
                        study_response=requests.get(url,headers=headers)
                except:
                        raise ValueError('ERROR REACHING %s' % (url))

                if study_response.status_code!=200:
                        raise ValueError('ERROR w/ %s : Code %s' % (url,study_response.status_code))
                        exit(1)

                if study_id in study_response.json()['organizations']:
                        return(str(study_response.json()['id']))

        raise ValueError('ERROR w/ %s : %s study was not found' % (url,study_id))

def check_analysis_types(file_manager_url,study_id,token):
        analysis_types=[]
        required_analysis_fields={}

        print("Retrieving analysis Types")
        headers={
                "Authorization" : "Bearer %s" % token
        }
        limit=20
        url="%s/schemas?hideSchema=true&limit=%s&offset=0&unrenderedOnly=false" % (file_manager_url,str(limit))

        try:
                response=requests.get(url,headers=headers)
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
                url="%s/schemas/%s?unrenderedOnly=false" % (file_manager_url,analysis_type)
                try:
                        response=requests.get(url,headers=headers)
                except:
                        raise ValueError('ERROR REACHING %s' % (url))

                if response.status_code!=200:
                        raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
                        exit(1)
                required_analysis_fields[analysis_type]={}
                required_analysis_fields[analysis_type]['fields']={}
                required_analysis_fields[analysis_type]['dataTypes']=[]

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

def ingest_file_data(
        file_metadata,
        analysis_metadata,
        workflow_metadata
):
        print("Ingesting file Data")
        data={}
        for file in [file_metadata,analysis_metadata,workflow_metadata]:
                if file:
                        data[file.split("/")[-1].removesuffix(".tsv").lower()]=pd.read_csv(file,sep='\t')

        ###Check mandatory columns

        for mandatory_col in ["analysisType","submitter_analysis_id"]:
                if mandatory_col not in data['analysis'].columns.values.tolist():
                        raise ValueError('Missing required column \'%s\' in file %s' % (mandatory_col,analysis_metadata))

        if workflow_metadata:
            for mandatory_col in ["submitter_workflow_id"]:
                    if mandatory_col not in data['workflow'].columns.values.tolist():
                            print(data['workflow'])
                            raise ValueError('Missing required column \'%s\' in workflow %s' % (mandatory_col,workflow_metadata))
        for mandatory_col in ["fileName","dataType"]:
                if mandatory_col not in data['files'].columns.values.tolist():
                        raise ValueError('Missing required column \'%s\' in files %s' % (mandatory_col,file_metadata))

        return(data)

def ingest_clinical_data(
        schema,
        sample_metadata,
        specimen_metadata,
        experiment_metadata,
        read_group_metadata
):
        print("Ingesting clinical Data")
        data={}
        for entity in [
                sample_metadata,
                specimen_metadata,
                experiment_metadata,
                read_group_metadata
                ]:
                if entity:
                        data[entity.split("/")[-1].removesuffix(".tsv").lower()]=pd.read_csv(entity,sep='\t')
        
        relational_mappings={}

        for key in ["experiment","read_group","sample","specimen"]:
                relational_mappings[key.lower()]={}
                relational_mappings[key]['primary']=[]
                relational_mappings[key]['foreign']=[]
                relational_mappings[key]['present']=True if key in data.keys() else False
                for entity in schema:
                        ####Check unique Id
                        if key.lower()==entity['name'].lower():
                                ####Check unique Id
                                for field in entity['fields']:
                                        if field.get('unique'):
                                                if field['unique']==True:
                                                        if relational_mappings[key]['present']:
                                                            if field['name']  not in data[key].columns.values.tolist():
                                                                    raise ValueError('Missing Unique identifier key %s in %s' % (field['name'] ,key))
                                                                    exit(1)
                                                            else:
                                                                    relational_mappings[key]['primary'].append(field['name'].lower())
                                                        else:
                                                            relational_mappings[key]['primary'].append(field['name'].lower())
                                ####Check foreignkeys
                                if entity.get("restrictions").get("foreignKey"):
                                        for foreignKey in entity.get("restrictions").get("foreignKey"):
                                                for mappings in foreignKey['mappings']:
                                                        if relational_mappings[key]['present']:
                                                            if mappings["local"] not in data[key].columns.values.tolist():
                                                                    raise ValueError('Missing foreign key %s in %s' % (mappings["local"],key))
                                                            else:
                                                                    relational_mappings[key]['foreign'].append(
                                                                            {
                                                                                    "entity":foreignKey['schema'].lower(),
                                                                                    "foreign":mappings['local'].lower()
                                                                            }
                                                                            )
                                                        else:
                                                            relational_mappings[key]['foreign'].append(
                                                                    {
                                                                            "entity":foreignKey['schema'].lower(),
                                                                            "foreign":mappings['local'].lower()
                                                                    }
                                                                    )                                                            

        return(data,relational_mappings)

def split_to_channel(file_data,clinical_data,relational_mapping,analysis_types,study_id,data_directory):
    print("MAPPING CLINICAL AND FILE DATA AT ANALYSIS LEVEL")
    ###Establish per analysis
    analyses={analysis:{} for analysis in file_data['analysis']['submitter_analysis_id'].unique().tolist()}
    for analysis in analyses.keys():
        analyses[analysis]['analysis']={
            "data":file_data['analysis'].query("submitter_analysis_id==@analysis").reset_index(drop=True).copy(),
            "status":True,
            "submission": True,
            "comments":[]
        }

        ###Checking for duplicates
        if len(analyses[analysis]['analysis']['data'])>1:
            analyses[analysis]['analysis']['status']=False
            analyses[analysis]['analysis']['comments'].append("Multiple conflict records with same submttier_analysis_id %s" % analysis)

        ###Checking for missing required values
        instance_analysis_type=analyses[analysis]['analysis']['data'].loc[0,"analysisType"]
        for col_group in analysis_types[instance_analysis_type]['fields']:
            if col_group not in ['analysisType','files','workflow','dataTypes']:
                for col in  analysis_types[instance_analysis_type]['fields'][col_group]:
                    for ind in analyses[analysis]['analysis']['data'].index.values.tolist():
                        if pd.isna(analyses[analysis]['analysis']['data'].loc[ind,col]):
                            analyses[analysis]['analysis']['status']=False
                            analyses[analysis]['analysis']['comments'].append("Analysis %s is missing required col %s" % (analysis,col))

        ###Check for study ID
        for ind in analyses[analysis]['analysis']['data'].index.values.tolist():
            if analyses[analysis]['analysis']['data'].loc[ind,"studyId"]!=study_id:
                analyses[analysis]['analysis']['status']=False
                analyses[analysis]['analysis']['comments'].append("Analysis %s is references studyId %s instead of %s" % (analysis,analyses[analysis]['analysis']['data'].loc[ind,"studyId"],study_id))

        ###Map workflow to analysis
        workflow=analyses[analysis]['analysis']['data'].loc[0,"submitter_workflow_id"]

        if pd.isna(workflow) or "workflow" not in file_data.keys():
            analyses[analysis]['workflow']={
                "data":pd.DataFrame(columns=["submitter_workflow_id"]),
                "status":False,
                "submission":"workflow" in file_data.keys(),
                "comments":[]
            }            
        else:
            analyses[analysis]['workflow']={
                "data":file_data['workflow'].query("submitter_workflow_id==@workflow").reset_index(drop=True).dropna().copy(),
                "status":True,
                "submission":"workflow" in file_data.keys(),
                "comments":[]
            }

        ###Check for workflow duplicates
        if len(analyses[analysis]['workflow']['data'])>1:
            analyses[analysis]['workflow']['status']=False
            analyses[analysis]['workflow']['comments'].append('Multiple records using the same submitter_workflow_id %s' % workflow)

        ###Check for workflows for applicable analysisTypes
        if "workflow" not in analysis_types[instance_analysis_type]['fields'].keys() and len(analyses[analysis]['workflow']['data'])!=0:
            analyses[analysis]['workflow']['status']=False
            analyses[analysis]['workflow']['comments'].append('Submitter_workflow_id %s is not required for analysisType %s' % (workflow,instance_analysis_type))
            analyses[analysis]['analysis']['status']=False
            analyses[analysis]['analysis']['comments'].append('Submitter_workflow_id %s is not required for analysisType %s' % (workflow,instance_analysis_type))
        elif "workflow" in analysis_types[instance_analysis_type]['fields'].keys() and len(analyses[analysis]['workflow']['data'])==0:
            analyses[analysis]['workflow']['status']=False
            analyses[analysis]['workflow']['comments'].append('Submitter_workflow_id %s is required for analysisType %s' % (workflow,instance_analysis_type))
            analyses[analysis]['analysis']['status']=False
            analyses[analysis]['analysis']['comments'].append('Submitter_workflow_id %s is required for analysisType %s' % (workflow,instance_analysis_type))

        ###Map files to analysis
        analyses[analysis]['files']={
            "data":file_data['files'].query("submitter_analysis_id==@analysis").reset_index(drop=True).dropna().copy(),
            "status":True,
            "submission":True,
            "comments":[]
        }

        ###Check if files exists
        for ind in analyses[analysis]['files']['data'].index.values.tolist():
            if not os.path.exists("%s/%s" % (data_directory,analyses[analysis]['files']['data'].loc[ind,"fileName"])):
                analyses[analysis]['files']['status']=False
                analyses[analysis]['files']['comments'].append("File %s/%s could not be found" % (data_directory,analyses[analysis]['files']['data'].loc[ind,"fileName"]))

        ###Checking for missing file required values
        instance_analysis_type=analyses[analysis]['analysis']['data'].loc[0,"analysisType"]
        for col_group in analysis_types[instance_analysis_type]['fields']:
            if col_group=='files':
                for col in  analysis_types[instance_analysis_type]['fields'][col_group]:
                    for ind in analyses[analysis]['files']['data'].index.values.tolist():
                            if pd.isna(analyses[analysis]['files']['data'].loc[ind,col]):
                                analyses[analysis]['files']['status']=False
                                analyses[analysis]['files']['comments'].append("Missing required col %s" % col)
            ###Check for datatypes
            if col_group=='dataTypes':
                for ind in analyses[analysis]['files']['data'].index.values.tolist():
                    if analyses[analysis]['files']['data'].loc[ind,"dataTypes"] not in analysis_types[instance_analysis_type]['fields']['dataTypes']:
                            analyses[analysis]['files']['status']=False
                            analyses[analysis]['files']['comments'].append("File %s with Datatype %s not part of approved files : %s" % (analyses[analysis]['files']['data'].loc[ind,"fileName"],analyses[analysis]['files']['data'].loc[ind,"dataTypes"]),",".join(analysis_types[instance_analysis_type]['fields']['dataTypes']))


            for col in ['fileName',"fileMd5sum"]:
                if len(analyses[analysis]['files']['data'])==0:
                    analyses[analysis]['files']['status']=False
                    analyses[analysis]['files']['comments'].append("No files provided for analysis %s" % (analysis))
                else:
                    if len(analyses[analysis]['files']['data'].query("%s=='%s'" % (col,analyses[analysis]['files']['data'].loc[ind,col])))>1:
                        analyses[analysis]['files']['status']=False
                        analyses[analysis]['files']['comments'].append("Multiple records found for col %s %s" % (col,analyses[analysis]['files']['data'].loc[ind,col]))

            ###Update relational mapping to add analysis. For now all analyses share the same mapping, may change in the future
            for dependency in analysis_types[instance_analysis_type]['externalValidations']:
                placeholder="https://submission.pcgl-dev.cumulus.genomeinformatics.org/validator/category/1/entity/experiment/exists?organization={studyId}&value={value}"
                dependency['url']=placeholder
                entity=re.findall(r'(?<=entity\/)([^\/]+)(?=\/exists)',dependency['url'])[0]
                submitter_id=dependency['jsonPath']

                relational_mapping['analysis']={
                    "primary":["submitter_analysis_id"],
                    "foreign":[{"entity":entity,"foreign":submitter_id}],
                    "present":True
                }
          
        
        for entity in ["analysis","experiment","read_group","sample","specimen"]:
            ###Add entity if not present:
            foreign_key=relational_mapping[entity]['foreign'][0]['foreign']
            foreign_entity=relational_mapping[entity]['foreign'][0]['entity']
            primary_key=relational_mapping[entity]['primary'][0]

            if foreign_entity=='participant':
                continue


            if entity not in analyses[analysis].keys():
                if  relational_mapping[entity]['present']:
                    search_value=analyses[analysis][foreign_entity]['data'].loc[:,foreign_key].values.tolist()[0]
                    analyses[analysis][entity]={
                        "data":clinical_data[entity].query("%s==@search_value" % foreign_key),
                        "status":True, ###True for now, will verify if values exist and reset status downstream
                        "submission":True,
                        "comments":[]
                    }
                    print(analysis,entity,"IS PRESENT")
                else:
                    # analyses[analysis][entity]={
                    #     "data":pd.DataFrame(),
                    #     "status":True, ###True for now, will verify if values exist and reset status downstream
                    #     "submission":False,
                    #     "comments":[]
                    # }
                    print(analysis,entity,"IS ABSENT")
            if relational_mapping[foreign_entity]['present'] and entity in analyses[analysis].keys():
                if len(analyses[analysis][entity]['data'].loc[:,foreign_key].values.tolist())>0:
                    print(analysis,entity,foreign_entity,"DATA IS PRESENT")
                    search_value=analyses[analysis][entity]['data'].loc[:,foreign_key].values.tolist()[0]
                    analyses[analysis][foreign_entity]={
                        "data":clinical_data[foreign_entity].query("%s==@search_value" % foreign_key),
                        "status":True, ###True for now, will verify if values exist and reset status downstream
                        "submission":True,
                        "comments":[]
                    }


            elif foreign_entity in analyses[analysis].keys():
                print(analysis,entity,foreign_entity,"DATA PREVIOUSLY ADDED")
                continue
            else:
                print(analysis,entity,foreign_entity,"DATA IS NOT PRESENT NEED TO INFER")

                if analyses[analysis][entity]['submission'] and len(analyses[analysis][entity]['data'].loc[:,foreign_key].values.tolist())>0:
                    print(analysis,entity,foreign_entity,"CAN INFER TO A LIMITED DEGREE")
                    search_value=analyses[analysis][entity]['data'].loc[:,foreign_key].values.tolist()[0]
                    analyses[analysis][foreign_entity]={
                        "data":pd.DataFrame(
                            [[search_value,np.nan]],
                            columns=[foreign_key,relational_mapping[foreign_entity]['foreign'][0]['foreign']]
                        ),
                        "status":False,
                        "submission":False,
                        "comments":["No data submitted"]
                    } 
                else:
                    print(analysis,entity,foreign_entity,"NOT ENOUGH TO INFER")
                    analyses[analysis][foreign_entity]={
                    "data":pd.DataFrame(columns=[foreign_key,relational_mapping[foreign_entity]['foreign'][0]['foreign']]),
                    "status":False,
                    "submission":False,
                    "comments":["No data submitted"]
                    }

        ###Check if file_r1 and file_r2 in read_group use correct files
        if "read_group" in analyses[analysis].keys():
                for ind in analyses[analysis]['read_group']['data'].index.values.tolist():
                        for file_read in ["file_r1","file_r2"]:
                                if file_read in analyses[analysis]['read_group']['data'].columns.values.tolist():
                                        if not pd.isna(analyses[analysis]['read_group']['data'].loc[ind,file_read]):
                                                if analyses[analysis]['read_group']['data'].loc[ind,file_read] not in analyses[analysis]['files']['data']['fileName'].values.tolist():
                                                        analyses[analysis]["read_group"]["status"]=False
                                                        analyses[analysis]["read_group"]["comments"].append("Read group %s references unlisted file %s " % (analysis,analyses[analysis]['read_group']['data'].loc[ind,file_read]))      

    return(analyses,relational_mapping)

def query_clinical_validator(clinical_url,study_id,category_id,token,entity,value):
        headers={
                "Authorization" : "Bearer %s" % token
        }


        #https://submission.pcgl-dev.cumulus.genomeinformatics.org/validator/category/1/entity/participant/exists?organization=EXAMPLE-CA&value=DONOR_01
        url="%s/validator/category/%s/entity/%s/exists?organization=%s&value=%s" % (clinical_url,category_id,entity.lower(),study_id,value)

        try:
                response=requests.get(url,headers=headers)
        except:
                raise ValueError('ERROR REACHING %s' % (url))

        if response.status_code!=200 and response.status_code!=404:
                raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
                exit(1)

        if response.json()['message']=='Record found':
            return(True)
        else:
            return(False)

def verify_submitted_data(clinical_url,study_id,category_id,token,entity,value,data):
        print("Verifying if submitted data is consistent")
        headers={
                "Authorization" : "Bearer %s" % token
        }
        comments=[]
        url="%s/data/category/%s/organization/%s/query?entityName=%s" % (clinical_url,category_id,study_id,entity)
        payload={
                "op": "and",
                "content": [
                        {
                                "op": "in",
                                "content": {
                                        "fieldName": unique_id_col,
                                        "value": [search_value]
                                        }
                                }
                        ]
                }

        try:
                response=requests.post(url,json=payload,headers=headers)
        except:
                raise ValueError('ERROR REACHING %s' % (url))

        if response.status_code!=200 and response.status_code!=404:
                raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
                exit(1)
        
        for ind in data.index.values.tolist():
                for col in data.columns.values.tolist():
                        valA=response.json()['records'][0][col]
                        valB=data.loc[ind,col]

                        if valA!=valB:
                                comments.append("Field '%s' is not consistent for record %s in entity %s. Specified - %s vs Comitted - %s" % (col,value,entity,valA,valB))
                                
        if len(comments)>0:
                return(True,comments)
        else:
                return(True,comments)

def check_registration(analysis_channels,clinical_url,category_id,study_id,token,relational_mappings):
    print("Verifying non-submitted IDs do not exist and dependencies are met")

    for analysis in analysis_channels:
        for entity in ["analysis","experiment","read_group","specimen","sample"]:
            if entity not in analysis_channels[analysis].keys():
                continue
            ##Check if already submitted
            if analysis_channels[analysis][entity]['submission']==True:
                primary_key=relational_mappings[entity]['primary'][0]
                for value in analysis_channels[analysis][entity]['data'][primary_key].values.tolist():
                ###Skip these as we only care about clinical

                    ###CHECK PRIMARY KEYS to make sure they arent submitted
                    ###Skip this check for analysis

                    if entity!='analysis':
                        verify_unsubmitted=query_clinical_validator(clinical_url,study_id,category_id,token,entity,value)

                        if verify_unsubmitted:
                                analysis_channels[analysis][entity]['comments'].append("The value '%s' of entity '%s' in study '%s' has already been registered" % (value,entity,study_id))
                                verification,comment=verify_submitted_data(
                                        clinical_url,
                                        study_id,
                                        category_id,
                                        token,
                                        entity,
                                        value,
                                        analysis_channels[analysis][entity]['data'].query("%s=='%s'" % (primary_key,value))
                                        )

                                if verification:
                                        analysis_channels[analysis][entity]['comments'].append(comment)
                    ###Check foreign keys to make sure dependency is met first by checking other submission
                    foreign_entity=relational_mappings[entity]['foreign'][0]['entity']
                    foreign_key=relational_mappings[entity]['foreign'][0]['foreign']

                    if foreign_entity in analysis_channels[analysis].keys():
                        for foreign_value in analysis_channels[analysis][entity]['data'][foreign_key]:
                            if foreign_value in analysis_channels[analysis][foreign_entity]['data'][foreign_key].unique().tolist():

                                if analysis_channels[analysis][foreign_entity]['status']:
                                    analysis_channels[analysis][entity]['comments'].append("Foreign value '%s' of entity '%s' in study '%s' will be registered" % (foreign_value,foreign_entity,study_id))
                                else:
                                    analysis_channels[analysis][entity]['status']=False
                                    analysis_channels[analysis][entity]['comments'].append("Foreign value '%s' of entity '%s' in study '%s' has not been registered" % (foreign_value,foreign_entity,study_id))                              
                            else:
                                ###Not found within submissions, check external DB
                                #verify_foreign_unsubmitted=query_clinical_validator(clinical_url,study_id,category_id,token,foreign_entity,val)
                                if query_clinical_validator(clinical_url,study_id,category_id,token,foreign_entity,foreign_value):
                                    analysis_channels[analysis][entity]['status']=True
                                    analysis_channels[analysis][entity]['comments'].append("Foreign value '%s' entity %s in study %s has been registered" % (foreign_value,foreign_entity,study_id))
                                else:
                                    analysis_channels[analysis][entity]['status']=False
                                    analysis_channels[analysis][entity]['comments'].append("Foreign value '%s' entity %s in study %s has not been registered" % (foreign_value,foreign_entity,study_id))
                    else:
                        for foreign_value in analysis_channels[analysis][entity]['data'][foreign_key]:
                            #verify_foreign_unsubmitted=query_clinical_validator(clinical_url,study_id,category_id,token,foreign_entity,val)
                            if query_clinical_validator(clinical_url,study_id,category_id,token,foreign_entity,foreign_value):
                                analysis_channels[analysis][entity]['comments'].append("Foreign value '%s' entity %s in study %s has been registered" % (foreign_value,foreign_entity,study_id))
                            else:
                                analysis_channels[analysis][entity]['status']=False
                                analysis_channels[analysis][entity]['comments'].append("Foreign value '%s' entity %s in study %s has not been registered" % (foreign_value,foreign_entity,study_id))                                                                                

        return(analysis_channels)

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

def save_outputs(analysis_channels,output_directory):

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    for analysis in analysis_channels.keys():
        if not os.path.isdir("%s/%s" % (output_directory,analysis)):
            os.makedirs("%s/%s" % (output_directory,analysis))

        check_flag=[]
        for entity in analysis_channels[analysis]:
            if analysis_channels[analysis][entity]['submission']:
                if len(analysis_channels[analysis][entity]['data'])>0:
                    if analysis_channels[analysis][entity]['status']:
                        tmp=analysis_channels[analysis][entity]['data'].copy()
                        tmp['comments']=";".join(list(set(analysis_channels[analysis][entity]['comments'])))
                        tmp.to_csv("%s/%s/registerable_%s.tsv" % (output_directory,analysis,entity),sep='\t',index=False)
                        check_flag.append(True)
                    else:
                        tmp=analysis_channels[analysis][entity]['data'].copy()
                        tmp['comments']=";".join(list(set(analysis_channels[analysis][entity]['comments'])))
                        tmp.to_csv("%s/%s/unregisterable_%s.tsv" % (output_directory,analysis,entity),sep='\t',index=False)
                        check_flag.append(False)

        if False in check_flag:
                file_name="%s/%s/%s.FAILURE" % (output_directory,analysis,analysis)
                with open(file_name, 'a') as file:
                        for entity in analysis_channels[analysis]:
                                file.write("%s\tstatus:%s\tsubmission:%s\n" % (entity,str(analysis_channels[analysis][entity]['status']),str(analysis_channels[analysis][entity]['submission'])))
                                for comment in analysis_channels[analysis][entity]['comments']:
                                        file.write(comment+"\n")
        else:
                file_name="%s/%s/%s.SUCCESS" % (output_directory,analysis,analysis)
                with open(file_name, 'a') as file:
                        for entity in analysis_channels[analysis]:
                                file.write("%s\tstatus:%s\tsubmission:%s\n" % (entity,str(analysis_channels[analysis][entity]['status']),str(analysis_channels[analysis][entity]['submission'])))
                                for comment in analysis_channels[analysis][entity]['comments']:
                                        file.write(comment+"\n")



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
        if args.output_directory: print("input:",args.output_directory)
        if args.data_directory: print("input:",args.data_directory)

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
        #check_clinical_study(
        #        args.clinical_url,
        #        args.study_id,
        #        args.token
        #        )

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

        file_data=ingest_file_data(
                args.file_metadata,
                args.analysis_metadata,
                args.workflow_metadata
        )       

        ###What do we have to register?
        ###R1d The pipeline shall query the clinical-submission service to check for the existence of the dependent entities for the submission. - Dependent entities registration is assumed pre-existing; failure here should skip the submission. The check should account for analysis specific dependencies.
        clinical_data,relational_mapping=ingest_clinical_data(
                clinical_schema,
                args.sample_metadata,
                args.specimen_metadata,
                args.experiment_metadata,
                args.read_group_metadata
                )
        
        analysis_channels,updated_relational_mapping=split_to_channel(
            file_data,
            clinical_data,
            relational_mapping,
            analysis_types,
            args.study_id,
            args.data_directory
        )

        registration_updated_analysis_channels=check_registration(
            analysis_channels,
            args.clinical_url,
            category_id,
            args.study_id,
            args.token,
            updated_relational_mapping
        )
        save_outputs(registration_updated_analysis_channels,args.output_directory)

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description='Tool: Check Clinical dependencies')
        parser.add_argument("-fi", "--file_metadata", dest="file_metadata", required=True, help="file metadata tsv")
        parser.add_argument("-an", "--analysis_metadata", dest="analysis_metadata", required=True, help="analysis metadata tsv")
        parser.add_argument("-wo", "--workflow_metadata", dest="workflow_metadata", required=False, help="workflow metadata tsv")
        parser.add_argument("-sa", "--sample_metadata", default=False, dest="sample_metadata", required=False, help="sample metadata tsv")
        parser.add_argument("-sp", "--specimen_metadata", default=False, dest="specimen_metadata", required=False, help="specimen metadata tsv")
        parser.add_argument("-ex", "--experiment_metadata", default=False, dest="experiment_metadata", required=False, help="experiment metadata tsv")
        parser.add_argument("-rg", "--read_group_metadata", default=False, dest="read_group_metadata", required=False, help="read_group metadata tsv")
        parser.add_argument("-cu", "--clinical_url", dest="clinical_url", required=True, help="Clinical URL")
        parser.add_argument("-fm", "--file_manager_url", dest="file_manager_url", required=True, help="File Manager URL")
        parser.add_argument("-si", "--study_id", dest="study_id", required=True, help="study_id")
        parser.add_argument("-t", "--token", dest="token", required=True, help="token")
        parser.add_argument("-od", "--output-directory", dest="output_directory", required=False, help="output directory where entity files are saved by analysis",default="output")
        parser.add_argument("-dd", "--data-directory", dest="data_directory", required=False, help="data directory where entity files are saved by analysis",default=os.getcwd())

        args = parser.parse_args()

        main(args)