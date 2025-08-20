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

def check_analysis_duplicates(analysis):
   if len(analysis['analysis'])>1:
      analysis['status']=False
      analysis['comments'].append("Multiple entries with the same submitter_analysis_id %s detected" % ";".join(analysis['analysis']['submitter_analysis_id'].unqiue().tolist()))
def check_workflows_duplicates(analysis):
   if "workflow" in analysis.keys():
      if len(analysis['workflow'])>1:
         analysis['status']=False
         analysis['comments'].append("Multiple entries with the same submitter_workflow_id %s detected" % ";".join(analysis['workflow']['submitter_workflow_id'].unqiue().tolist()))
def check_file_filename_duplicates(analysis):
   if len(analysis['files'].groupby("fileName").count().query('fileMd5sum>1'))>1:
      for fileName in analysis['files'].groupby("fileName").count().query('fileMd5sum>1').index.values.tolist():
         analysis['status']=False
         analysis['comments'].append("Multiple entries with the same fileName %s detected" % fileName)
def check_file_filename_filem5d(analysis):
   if len(analysis['files'].groupby("fileMd5sum").count().query('fileName>1'))>1:
      for fileName in analysis['files'].groupby("fileMd5sum").count().query('fileName>1').index.values.tolist():
         analysis['status']=False
         analysis['comments'].append("Multiple entries with the same fileMd5sum %s detected" % fileMd5sum)
def check_study_id(analysis,study_id):
   for ind in analysis['analysis'].index.values.tolist():
      if analysis['analysis'].loc[ind,"studyId"] != study_id:
         analysis['status']=False
         analysis['comments'].append("Invalid studyId record found for analysis %s. Expected %s found %s" % (analysis['analysis'].loc[ind,"studyId"],study_id,analysis['analysis'].loc[ind,"studyId"]))       
def check_file_minimum(analysis):
   if len(analysis['files'])==0:
      analysis['status']=False
      analysis['comments'].append("At minimum 1 file record expected")       
def check_file_exists(analysis,data_directory):
   for ind in analysis['files'].index.values.tolist():
      if not os.path.exists("%s/%s" % (data_directory,analysis['files'].loc[ind,"fileName"])):
            analysis['status']=False
            analysis['comments'].append("File %s/%s could not be found" % (data_directory,analysis['files'].loc[ind,"fileName"]))
def check_workflow_analysis(analysis,analysis_types):
   analysis_types_w_workflows=[]

   for schema in analysis_types:
      if "workflow" in analysis_types[schema]['fields']:
         analysis_types_w_workflows.append(schema)

   if analysis.get('analysis')["analysisType"].values.tolist()[0] in analysis_types_w_workflows and not analysis.get('workflows'):
      analysis['status']=False
      analysis['comments'].append("AnalysisType %s expects workflow info" % (analysis.get('analysis')["analysisType"].values.tolist()[0]))
   elif analysis.get('analysis')["analysisType"].values.tolist()[0] not in analysis_types_w_workflows and analysis.get('workflows'):
      analysis['status']=False
      analysis['comments'].append("AnalysisType %s does not require workflow info, but workflow record found" % (analysis.get('analysis')["analysisType"].values.tolist()[0]))      

def check_workflow_datatypes(analysis,analysis_types):
   analysisType=analysis.get('analysis')["analysisType"].values.tolist()[0]

   for ind in analysis.get('files').index.values.tolist():
      if analysis.get('files').loc[ind,'fileType'] not in analysis_types[analysisType]['dataTypes']:
         analysis['status']=False
         analysis['comments'].append("File %s of fileType %s is not accepted for analysisType %s. Allowed fileTypes are : %s" %  analysis.get('files').loc[ind,'fileName'],analysis.get('files').loc[ind,'fileType'],analysisType,",".join(analysis_types[analysisType]['dataTypes']))
def check_read_group_exists(analysis,analysisTypes,token,clinical_url,category_id,study_id):
   analysisType=analysis.get('analysis')["analysisType"].values.tolist()[0]
   experiment_id=analysis.get('analysis')["submitter_experiment_id"].values.tolist()[0]
   internal_check=False
   external_check=False
   if analysisType!='sequenceExperiment' and analysis.get('read_groups'):
      analysis['status']=False
      analysis['comments'].append("AnalysisType %s does not require read_group info, but read_group records found" % (analysis.get('analysis')["analysisType"].values.tolist()[0]))
   
   if analysisType=='sequenceExperiment':
      if analysis.get('read_groups'):
         internal_check=True
      
      headers={
               "Authorization" : "Bearer %s" % token
      }
      comments=[]
      url="%s/data/category/%s/organization/%s/query?entityName=%s" % (clinical_url,category_id,study_id,"read_group")
      payload={
               "op": "and",
               "content": [
                     {
                              "op": "in",
                              "content": {
                                       "fieldName": "submitter_experiment_id",
                                       "value": [experiment_id]
                                       }
                              }
                     ]
               }

      try:
         response=requests.post(url,json=payload,headers=headers)
      except:
         analysis['status']=False
         analysis['comments'].append("check_read_group_exists - ERROR REACHING %s" % (url))

      if response.status_code!=200 and response.status_code!=404:
         analysis['status']=False
         analysis['comments'].append("check_read_group_exists - ERROR w/ %s : Code %s" % (url,response.status_code))
         return
      elif response.status_code==200:
         external_check=True
      elif response.status_code==404:
         external_check=False
      

      if not internal_check and not external_check:
         analysis['status']=False
         analysis['comments'].append("check_read_group_exists - No read_groups related to submitter_experiment_id %s found for analysisType %s" % (experiment_id,analysisType))

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

        if response.json()['message']=='Record found':
            return(True)
        else:
            return(False)
def query_registered_data(token,clinical_url,category_id,entity,study_id,primary_key,ind,analysis):
   print("Verifying if submitted data is consistent")
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
                                    "value": [analysis[entity].loc[ind,primary_key]]
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

   for col in analysis[entity].columns.values.tolist():
      if not pd.isna(analysis[entity].loc[ind,col]):
         if response.json()['records'][0]['data'].get(col):
            valA=response.json()['records'][0]['data'][col]

            valB=analysis[entity].loc[ind,col]

            if valA!=valB:
               analysis["comments"].append("Field '%s' is not consistent for record %s in entity %s. Specified - %s vs Comitted - %s" % (col,analysis[entity].loc[ind,primary_key],entity,valA,valB))
               analysis['status']=False
         else:
            analysis["comments"].append("New Field '%s' detected for existing record %s in entity %s. Contact PCGL Admin to perform this update seperately" % (col,analysis[entity].loc[ind,primary_key],entity))
            analysis['status']=False

def check_registered_entities(analysis,clinical_url,category_id,study_id,relational_mapping,token):
   for entity in ["read_group","specimen","sample","experiment"]:
      if entity in analysis:
         for primary_key in relational_mapping[entity]['primary']:
            for ind in analysis[entity].index.values.tolist():
               url="%s/validator/category/%s/entity/%s/exists?organization=%s&value=%s" % (clinical_url,category_id,entity,study_id,analysis[entity].loc[ind,primary_key])

               if query_clinical_validator(url,token):
                  print(url)
                  query_registered_data(token,clinical_url,category_id,entity,study_id,primary_key,ind,analysis)


def main(args):
   if args.file_metadata: print("input:",args.file_metadata)
   if args.analysis_metadata: print("input:",args.analysis_metadata)
   if args.workflow_metadata: print("input:",args.workflow_metadata)
   if args.sample_metadata: print("input:",args.sample_metadata)
   if args.specimen_metadata: print("input:",args.specimen_metadata)
   if args.experiment_metadata: print("input:",args.experiment_metadata)
   if args.read_group_metadata: print("input:",args.read_group_metadata)
   if args.relational_mapping: print("input:",args.relational_mapping)
   if args.analysis_types: print("input:",args.analysis_types)
   if args.clinical_url: print("input:",args.clinical_url)
   if args.file_manager_url: print("input:",args.file_manager_url)
   if args.study_id: print("input:",args.study_id)
   if args.token: print("input:",args.token)
   if args.data_directory: print("input:",args.data_directory)

   with open(args.relational_mapping, 'r') as file:
      relational_mapping = json.load(file)

   with open(args.analysis_types, 'r') as file:
      analysis_types = json.load(file)

   category_id=retrieve_category_id(
      args.clinical_url,
      args.study_id,
      args.token
   )

   analysis={
      "status":True,
      "comments":[]
   }

   for metadata,key in zip(
        [args.analysis_metadata,args.file_metadata,args.workflow_metadata,args.sample_metadata,args.specimen_metadata,args.experiment_metadata,args.read_group_metadata],
        ["analysis","files","workflow","sample","specimen","experiment","read_group"],
   ):
      if metadata:
         analysis[key]=pd.read_csv(metadata,sep='\t')

   ###Flag duplicate analyses
   check_analysis_duplicates(analysis)

   ###Check for workflow duplicates
   check_workflows_duplicates(analysis)

   ###Check for duplicate ['fileName',"fileMd5sum"]:
   check_file_filename_duplicates(analysis)
   check_file_filename_filem5d(analysis)
   ###Check for study ID
   check_study_id(analysis,args.study_id)
   ###Check if files exists
   check_file_exists(analysis,args.data_directory)
   ###Check workflows map to correct analysisType
   check_workflow_analysis(analysis,analysis_types)
   ###Check for datatypes
   check_workflow_datatypes(analysis,analysis_types)

   #check_read_group_exists(analysis,analysis_types,args.token,args.clinical_url,category_id,args.study_id)
   ###If entity pre-registered check values
   check_registered_entities(analysis,args.clinical_url,category_id,args.study_id,relational_mapping,args.token)

   if len(analysis['comments'])>1:
      raise ValueError(".\n".join(analysis['comments']))


   

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description='Tool: Check Clinical dependencies')

   parser.add_argument("-fi", "--file_metadata", dest="file_metadata", required=True, help="file metadata tsv")
   parser.add_argument("-an", "--analysis_metadata", dest="analysis_metadata", required=True, help="analysis metadata tsv")
   parser.add_argument("-wo", "--workflow_metadata", default=False,dest="workflow_metadata", required=False, help="workflow metadata tsv")
   parser.add_argument("-sa", "--sample_metadata", default=False, dest="sample_metadata", required=False, help="sample metadata tsv")
   parser.add_argument("-sp", "--specimen_metadata", default=False, dest="specimen_metadata", required=False, help="specimen metadata tsv")
   parser.add_argument("-ex", "--experiment_metadata", default=False, dest="experiment_metadata", required=False, help="experiment metadata tsv")
   parser.add_argument("-rg", "--read_group_metadata", default=False, dest="read_group_metadata", required=False, help="read_group metadata tsv")
   parser.add_argument("-rm", "--relational_mapping", dest="relational_mapping", required=True, help="relational mapping json")
   parser.add_argument("-at", "--analysis_types", dest="analysis_types", required=True, help="analysis types json")
   parser.add_argument("-cu", "--clinical_url", dest="clinical_url", required=True, help="Clinical URL")
   parser.add_argument("-fm", "--file_manager_url", dest="file_manager_url", required=True, help="File Manager URL")
   parser.add_argument("-si", "--study_id", dest="study_id", required=True, help="study_id")
   parser.add_argument("-t", "--token", dest="token", required=True, help="token")
   parser.add_argument("-od", "--output_directory", default="output", dest="output_directory", required=False, help="output directory")
   parser.add_argument("-dd", "--data_directory", dest="data_directory", required=True, help="data directory where entity files are saved by analysis",default=os.getcwd())
   args = parser.parse_args()

   main(args)