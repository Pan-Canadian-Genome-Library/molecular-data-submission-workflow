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
                        #print(analysis_type,required,"string" if response.json()['schema']['properties'][required]=='string' else "?????")
                        #if response.json()['schema']['properties']
                        #required_analysis_fields[analysis_type][required]=[]
                
        return(required_analysis_fields)
                
def check_analysis(analysis_types,analysis_metadata_csv):
        print(analysis_types)


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

        for key in [key.lower() for key in data.keys()]:
                relational_mappings[key.lower()]={}
                relational_mappings[key]['primary']=[]
                relational_mappings[key]['foreign']=[]
                for entity in schema:
                        ####Check unique Id
                        if key.lower()==entity['name'].lower():
                                ####Check unique Id
                                for field in entity['fields']:
                                        if field.get('unique'):
                                                if field['unique']==True:
                                                        if field['name']  not in data[key].columns.values.tolist():
                                                                raise ValueError('Missing Unique identifier key %s in %s' % (field['name'] ,key))
                                                                exit(1)
                                                        else:
                                                                #print("example",key,field['name'],data[key][field['name']].values.tolist())
                                                                relational_mappings[key]['primary'].append(field['name'].lower())
                                ####Check foreignkeys
                                if entity.get("restrictions").get("foreignKey"):
                                        for foreignKey in entity.get("restrictions").get("foreignKey"):
                                                for mappings in foreignKey['mappings']:
                                                        if mappings["local"] not in data[key].columns.values.tolist():
                                                                raise ValueError('Missing foreign key %s in %s' % (mappings["local"],key))
                                                        else:
                                                                relational_mappings[key]['foreign'].append(
                                                                        {
                                                                                "entity":foreignKey['schema'].lower(),
                                                                                "foreign":mappings['local'].lower()
                                                                        }

                                                                        )

        return(data,relational_mappings)

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
                        #print(file)
                        #print(data[file.split("/")[-1].removesuffix(".tsv").lower()])
        ###Check mandatory columns
        #print(data['analysis'].columns.values.tolist())
        for mandatory_col in ["analysisType","submitter_analysis_id"]:
                if "analysisType" not in data['analysis'].columns.values.tolist():
                        raise ValueError('Missing required column \'%s\' in file %s' % (mandatory_col,analysis_metadata))

        return(data)




def check_clinical_registration(clinical_url,data,relational_mapping,category_id,study_id,token):
        print("Check current data and seperated based on registration")
        registered_data={}
        unregistered_data={}

        headers={
                "Authorization" : "Bearer %s" % token
        }
        ###Copy over dataframes
        for entity in data.keys():
                registered_data[entity]=pd.DataFrame(columns=data[entity].columns.values.tolist())
                unregistered_data[entity]=pd.DataFrame(columns=data[entity].columns.values.tolist())


        for entity in relational_mapping.keys():

                for ind in data[entity].index.values.tolist():
                ###Check if registered
                        for primary in relational_mapping[entity]['primary']:
                                #https://submission.pcgl-dev.cumulus.genomeinformatics.org/validator/category/1/entity/participant/exists?organization=EXAMPLE-CA&value=DONOR_01
                                url="%s/validator/category/%s/entity/%s/exists?organization=%s&value=%s" % (clinical_url,category_id,entity.lower(),study_id,data[entity].loc[ind,primary])

                                try:
                                        response=requests.get(url,headers=headers)
                                except:
                                        raise ValueError('ERROR REACHING %s' % (url))

                                if response.status_code!=200 and response.status_code!=404:
                                        raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
                                        exit(1)

                                if response.json()['message']=='Record found':
                                        registered_data[entity].loc[len(registered_data[entity])]=data[entity].loc[ind,:].values.tolist()
                                else:
                                        unregistered_data[entity].loc[len(unregistered_data[entity])]=data[entity].loc[ind,:].values.tolist()
        

                                

        return(registered_data,unregistered_data)

def check_unregistered_clinical(schema,clinical_url,unregistered,relational_mapping,category_id,study_id,token):
        print("Check unregistered data's dependency")
        entity_order=["specimen","sample","experiment","read_group"]
        unregisterable_data={}
        registerable_data={}

        headers={
                "Authorization" : "Bearer %s" % token
        }
        ###Copy over dataframes
        for entity in unregistered.keys():
                unregisterable_data[entity]=pd.DataFrame(columns=unregistered[entity].columns.values.tolist())
                unregisterable_data[entity]["reason"]=None
                registerable_data[entity]=pd.DataFrame(columns=unregistered[entity].columns.values.tolist())

        tmp=entity_order.copy()
        for entity in tmp:
                if entity not in unregistered.keys():
                        entity_order.remove(entity)
        
        for entity in entity_order:

                mandatory_fields=[]
                for lectern_entity in schema:
                        if entity.lower()==lectern_entity['name']:
                                for field in lectern_entity['fields']:
                                        if field.get('restrictions'):
                                                if field.get('restrictions').get('required'):
                                                        if field['restrictions']['required']:
                                                                mandatory_fields.append(field['name'])

                for ind in unregistered[entity].index.values.tolist():
                        print(unregistered)
                        print("AAAAAA",entity)

                        for field in mandatory_fields:
                                if unregistered[entity].loc[ind,field]==None or unregistered[entity].loc[ind,field]==np.NaN or pd.isna(unregistered[entity].loc[ind,field]):
                                        unregisterable_data[entity].loc[len(unregisterable_data[entity])]=unregistered[entity].loc[ind,:].values.tolist()+["Missing required field %s" % field]



                        for foreign in relational_mapping[entity]['foreign']:
                                ###Check if dependency is registered
                                url="%s/validator/category/%s/entity/%s/exists?organization=%s&value=%s" % (clinical_url,category_id,foreign['entity'],study_id,unregistered[entity].loc[ind,foreign['foreign']])


                                try:
                                        response=requests.get(url,headers=headers)
                                except:
                                        raise ValueError('ERROR REACHING %s' % (url))

                                if response.status_code!=200 and response.status_code!=404:
                                        raise ValueError('ERROR w/ %s : Code %s' % (url,response.status_code))
                                        exit(1)

                                print("WTF",unregistered)
                                if response.json()['message']=='Record found':
                                        registerable_data[entity].loc[len(registerable_data[entity])]=unregistered[entity].loc[ind,:].values.tolist()
                                ###If not registered check if dependency is in the to register list
                                elif unregistered.get(foreign['entity']):
                                        if unregistered[entity].loc[ind,foreign['foreign']] in unregistered[foreign['entity']][foreign['foreign']].values.tolist():
                                                registerable_data[entity].loc[len(registerable_data[entity])]=unregistered[entity].loc[ind,:].values.tolist()
                                        else:
                                                unregisterable_data[entity].loc[len(unregisterable_data[entity])]=unregistered[entity].loc[ind,:].values.tolist()+["Missing dependency %s" % foreign['foreign']]              
                                else:
                                        unregisterable_data[entity].loc[len(unregisterable_data[entity])]=unregistered[entity].loc[ind,:].values.tolist()+["Missing dependency %s" % foreign['foreign']]
                
        return(unregisterable_data,registerable_data)

                                
def check_registered_clinical(clinical_url,registered,relational_mapping,category_id,study_id,token):
        headers={
                "Authorization" : "Bearer %s" % token
        }
        print("Check registered data is consistent")
        consistent_registered_data={}
        inconsistent_registered_data={}
        for entity in registered.keys():
                inconsistent_registered_data[entity]=pd.DataFrame(columns=registered[entity].columns.values.tolist())
                inconsistent_registered_data[entity]["reason"]=None
                consistent_registered_data[entity]=pd.DataFrame(columns=registered[entity].columns.values.tolist())

        return(consistent_registered_data,inconsistent_registered_data)
        for entity in registered.keys():
                #print(entity)
                #print(registered[entity].head())
                for ind in registered[entity].index.values.tolist():
                        unique_id_col=relational_mapping[entity]['primary'][0]
                        search_value=registered[entity].loc[ind,unique_id_col]

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

                        for col in registered[entity].loc[ind,:].dropna().index.tolist():
                                print(response.text)
                                valA=response.json()['records'][0][col]
                                valB=registered[entity].loc[ind,col]
                                if valA!=valB:
                                        unregisterable_data[entity].loc[len(unregisterable_data[entity])]=unregistered[entity].loc[ind,:].values.tolist()+["Missing required field %s" % field]
                                        inconsistent_registered_data[entity].loc[len(inconsistent_registered_data[entity])]=registered[entity].loc[ind,:].values.tolist()+["Field %s is not consistent. Specified - %s vs Comitted - %s" % (col,valA,valB) ]
                                        continue
                        
                        consistent_registered_data[entity].loc[len(consistent_registered_data[entity])]=registered[entity].loc[ind,:].values.tolist()

        return(consistent_registered_data,inconsistent_registered_data)

def check_registerable_analyses(file_data,category_id,study_id,token,inconsistent_registered_data,unregisterable_clinical_data,registerable_clinical_data):

        for entity in file_data.keys():
                print(entity)
                print(file_data[entity])
        #Per analysi, check for duplicate submitterIdanalysis
        #PEr analysis check if used before / No present
        #Per analysis check if required columns are present
        #Per Analysis check if files are present and if files are correct types
        #Per analys check if required workflow info is present
        #
def main(args):
        print("input:",args.sample_metadata)
        print("input:",args.specimen_metadata)
        print("input:",args.experiment_metadata)
        print("input:",args.read_group_metadata)
        print("input:",args.clinical_url)
        print("input:",args.study_id)
        print("input:",args.token)

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

        ###R1a - The pipeline shall verify the presence of required input metadata files and molecular files. Fundamental to ensure all inputs are available before proceeding.

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

        registered,unregistered=check_clinical_registration(
                args.clinical_url,
                clinical_data,
                relational_mapping,
                category_id,
                args.study_id,
                args.token)

        ###R1c - The pipeline shall query the clinical-submission service to check for registered participants. Participant registration is assumed pre-existing; failure here should halt the pipeline.
        unregisterable_clinical_data,registerable_clinical_data=check_unregistered_clinical(
                clinical_schema,
                args.clinical_url,
                unregistered,
                relational_mapping,
                category_id,
                args.study_id,
                args.token)

        consistent_registered_data,inconsistent_registered_data=check_registered_clinical(
                args.clinical_url,
                registered,
                relational_mapping,
                category_id,
                args.study_id,
                args.token)

        check_registerable_analyses(
                file_data,
                category_id,
                args.study_id,
                args.token,
                inconsistent_registered_data,
                unregisterable_clinical_data,
                registerable_clinical_data
         )

        for entity in registerable_clinical_data.keys():
                if len(registerable_clinical_data[entity])>0:
                        registerable_clinical_data[entity].to_csv("registerable_%s.tsv" % (entity),sep='\t',index=False)
        for entity in unregisterable_clinical_data.keys():
                if len(unregisterable_clinical_data[entity])>0:
                        unregisterable_clinical_data[entity].to_csv("unregisterable_%s.tsv" % (entity),sep='\t',index=False)
        for entity in registered.keys():
                if len(registered[entity])>0:
                        registered[entity].to_csv("registered_%s.tsv" % (entity),sep='\t',index=False)




if __name__ == "__main__":
        parser = argparse.ArgumentParser(description='Tool: Check Clinical dependencies')
        parser.add_argument("-fi", "--file_metadata_csv", dest="file_metadata", required=True, help="file metadata csv")
        parser.add_argument("-an", "--analysis_metadata_csv", dest="analysis_metadata", required=True, help="analysis metadata csv")
        parser.add_argument("-wo", "--workflow_metadata_csv", dest="workflow_metadata", required=True, help="workflow metadata csv")
        parser.add_argument("-sa", "--sample_metadata_csv", default=False, dest="sample_metadata", required=False, help="sample metadata csv")
        parser.add_argument("-sp", "--specimen_metadata_csv", default=False, dest="specimen_metadata", required=False, help="specimen metadata csv")
        parser.add_argument("-ex", "--experiment_metadata_csv", default=False, dest="experiment_metadata", required=False, help="experiment metadata csv")
        parser.add_argument("-rg", "--read_group_metadata_csv", default=False, dest="read_group_metadata", required=False, help="read_group metadata csv")
        parser.add_argument("-cu", "--clinical_url", dest="clinical_url", required=True, help="Clinical URL")
        parser.add_argument("-fm", "--file_manager_url", dest="file_manager_url", required=True, help="File Manager URL")
        parser.add_argument("-si", "--study_id", dest="study_id", required=True, help="study_id")
        parser.add_argument("-t", "--token", dest="token", required=True, help="token")

        args = parser.parse_args()

        main(args)
