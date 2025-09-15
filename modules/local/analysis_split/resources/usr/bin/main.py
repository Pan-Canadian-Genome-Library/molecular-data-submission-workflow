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
import datetime

def split_analyses(analysis_df):
    print("Creating dictionary based off of available analyses")
    return(
        {
            analysis:{
                "analysis":{
                    'data':analysis_df.get('data').query("submitter_analysis_id==@analysis"),
                    "submitted":True
                    },
                "comments":[],
                "status":True} for analysis in analysis_df.get('data')['submitter_analysis_id'].unique().tolist()
        }
    )

def flag_duplicate_analyses(analyses):
    print("Checking for duplicate analyses")
    for analysis in analyses:
        if len(analyses.get(analysis).get("analysis").get('data'))>1:
            analyses[analysis]['status']=False
            analyses[analysis]['comments'].append("Duplicate analyses found for %s : %s conflicting records" % (analysis,str(len(analyses.get(analysis).get("analysis")))))

def map_analysis_dependencies(analyses,relational_mapping,data):
    print("Updating relational_mapping to include analysis foreign keys")
    for analysis in analyses:
        for analysisType in analyses.get(analysis).get("analysis").get('data')['analysisType'].unique().tolist():
            foreign_entity=relational_mapping.get("analysis").get("analysisTypes").get(analysisType).get("foreign").get("entity")
            foreign_key=relational_mapping.get("analysis").get("analysisTypes").get(analysisType).get("foreign").get("foreign")
            foreign_values=analyses.get(analysis).get('analysis').get('data').loc[:,foreign_key].values.tolist()

            analyses[analysis][foreign_entity]={}
            if foreign_entity in data.keys() and foreign_entity!='participant':
                analyses[analysis][foreign_entity]['data']=data[foreign_entity].get('data').query("%s==@foreign_values" % foreign_key)
                analyses[analysis][foreign_entity]['submitted']=True
            else:
                analyses[analysis][foreign_entity]['data']=pd.DataFrame(foreign_values if isinstance(foreign_values, list) else [foreign_values],columns=[foreign_key])
                analyses[analysis][foreign_entity]['submitted']=False

def map_biospecimen_entities(analyses,relational_mapping,data):
    print("Mapping biospecimens to analyses")
    for entity in ["analysis","experiment","sample","specimen","read_group"]:
        for analysis in analyses:
            if entity=='analysis':
                analysis_type=analyses.get(analysis).get(entity).get('data').loc[:,"analysisType"].values.tolist()[0]
                #print(entity,analysis_type)
                foreign_entity=relational_mapping.get(entity).get("analysisTypes").get(analysis_type).get("foreign").get("entity")
                foreign_key=relational_mapping.get(entity).get("analysisTypes").get(analysis_type).get("foreign").get("foreign")
            else:
                foreign_entity=relational_mapping.get(entity).get("foreign")[0].get("entity")
                foreign_key=relational_mapping.get(entity).get("foreign")[0].get("foreign")
            ###Only check if biospecimen records are local
            #print(analyses.get(analysis))
            if analyses.get(analysis).get(entity):
                #If experiment exists, use experiment to find sample
                if analyses.get(analysis).get(entity).get('submitted') and foreign_entity!='participant':
                    #print(entity,foreign_entity)
                    foreign_values=analyses.get(analysis).get(entity).get('data').loc[:,foreign_key].values.tolist()
                    analyses[analysis][foreign_entity]={}
                    if foreign_entity in data.keys():
                        analyses[analysis][foreign_entity]['data']=data[foreign_entity].get('data').query("%s==@foreign_values" % foreign_key)
                        analyses[analysis][foreign_entity]['submitted']=data.get(foreign_entity).get('submitted')
                    else:
                        analyses[analysis][foreign_entity]['data']=pd.DataFrame(foreign_values if isinstance(foreign_values, list) else [foreign_values],columns=[foreign_key])
                        analyses[analysis][foreign_entity]['submitted']=False
                # else:
                #     analyses[analysis][foreign_entity]['data']=pd.DataFrame([foreign_values],columns=[foreign_key])
                #     analyses[analysis][foreign_entity]['submitted']=False                   
            elif analyses.get(analysis).get(foreign_entity):

                #From the other end, if experiment exists, use experiment to find read_group
                if analyses.get(analysis).get(foreign_entity).get('submitted') and foreign_entity!='participant':
                    foreign_values=analyses.get(analysis).get(foreign_entity).get('data').loc[:,foreign_key].values.tolist()
                    analyses[analysis][entity]={}
                    if data.get(entity):
                        if data.get('submitted'):
                            analyses[analysis][entity]['data']=data[entity].get('data').query("%s==@foreign_values" % foreign_key)
                            analyses[analysis][entity]['submitted']=True
                        else:
                            analyses[analysis][entity]['data']=pd.DataFrame(foreign_values if isinstance(foreign_values, list) else [foreign_values],columns=[foreign_key])
                            analyses[analysis][entity]['submitted']=False
                    else:
                        analyses[analysis][entity]['data']=pd.DataFrame(foreign_values if isinstance(foreign_values, list) else [foreign_values],columns=[foreign_key])
                        analyses[analysis][entity]['submitted']=False
            else:
                ###
                print("Nothing to map for entity %s in analysis %s" % (entity,analysis))

def map_files(analyses,relational_mapping,data):
    print("Mapping files to analyses")
    entity="files"
    foreign_entity=relational_mapping.get(entity).get("foreign")[0].get("entity")
    foreign_key=relational_mapping.get(entity).get("foreign")[0].get("foreign")
    for analysis in analyses:
        analyses[analysis][entity]={}
        foreign_values=analyses.get(analysis).get(foreign_entity).get('data').loc[:,foreign_key].values.tolist()
        analyses[analysis][entity]['data']=data[entity].get('data').query("%s==@foreign_values" % foreign_key)
        analyses[analysis][entity]['submitted']=True

def map_workflows(analyses,relational_mapping,data):
    print("Mapping workflows to analyses")
    entity="workflow"
    for analysis in analyses:
        analysis_types=analyses.get(analysis).get("analysis").get("data").loc[:,"analysisType"]
        for analysis_type in analysis_types:
            if analysis_type in relational_mapping.get(entity).get('analysisTypes').keys():
                analyses[analysis][entity]={}
                if entity in data:
                    analyses[analysis][entity]['data']=data[entity].get('data').query("submitter_analysis_id==@analysis")
                    analyses[analysis][entity]['submitted']=True

def save_outputs(analyses,output_directory):       
    print("Saving analysis specific TSVs and status.yml")
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    for analysis in analyses.keys():
        if not os.path.isdir("%s/%s" % (output_directory,analysis)):
            os.makedirs("%s/%s" % (output_directory,analysis))

        for entity in analyses[analysis].keys():
            if entity not in ['submitted','comments',"status"]:
                if analyses.get(analysis).get(entity).get('submitted'):
                    if len(analyses.get(analysis).get(entity).get('data'))>0:
                        analyses.get(analysis).get(entity).get('data').to_csv("%s/%s/%s.tsv" % (output_directory,analysis,entity),sep='\t',index=False)
        with open("%s/%s/%s_check_submission_dependencies_analysis_split_status.yml" % (output_directory,analysis,analysis), 'w') as file:
            file.write(
            """
process: "CHECK_SUBMISSION_DEPENDENCIES:ANALYSIS_SPLIT"
status: "%s"
exit_code: %s
timestamp: "%s"
work_directory: "%s"
details:
    analysis_id: "%s"
    error_message: "%s"

            """ % (
                    "PASS" if analyses.get(analysis).get("status") else "FAILED",
                    "0" if analyses.get(analysis).get("status") else "1",
                    datetime.datetime.now(),
                    os.getcwd(),
                    analysis,
                    "\n".join(analyses.get(analysis).get("comments") )
                )
            )



def main(args):
    if args.file_metadata: print("input:",args.file_metadata)
    if args.analysis_metadata: print("input:",args.analysis_metadata)
    if args.workflow_metadata: print("input:",args.workflow_metadata)
    if args.sample_metadata: print("input:",args.sample_metadata)
    if args.specimen_metadata: print("input:",args.specimen_metadata)
    if args.experiment_metadata: print("input:",args.experiment_metadata)
    if args.read_group_metadata: print("input:",args.read_group_metadata)
    if args.relational_mapping: print("input:",args.relational_mapping)

    with open(args.relational_mapping, 'r') as file:
            relational_mapping = json.load(file)

    data={}
    for metadata,key in zip(
        [args.analysis_metadata,args.file_metadata,args.workflow_metadata,args.sample_metadata,args.specimen_metadata,args.experiment_metadata,args.read_group_metadata],
        ["analysis","files","workflow","sample","specimen","experiment","read_group"],
    ):
        if metadata:
            data[key]={}
            data[key]['data']=pd.read_csv(metadata,sep='\t')
            data[key]['submitted']=True

    analyses=split_analyses(data['analysis'])

    flag_duplicate_analyses(analyses)

    map_analysis_dependencies(analyses,relational_mapping,data)

    ###This is currently making the assumption that all analyses are mapped to experiment, address this when this changes in the future
    map_biospecimen_entities(analyses,relational_mapping,data)

    map_files(analyses,relational_mapping,data)

    map_workflows(analyses,relational_mapping,data)
    save_outputs(analyses,args.output_directory)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Tool: Check Clinical dependencies')
    parser.add_argument("-fi", "--file_metadata", dest="file_metadata", required=True, help="file metadata tsv")
    parser.add_argument("-an", "--analysis_metadata", dest="analysis_metadata", required=True, help="analysis metadata tsv")
    parser.add_argument("-wo", "--workflow_metadata", dest="workflow_metadata", required=False, help="workflow metadata tsv")
    parser.add_argument("-sa", "--sample_metadata", default=False, dest="sample_metadata", required=False, help="sample metadata tsv")
    parser.add_argument("-sp", "--specimen_metadata", default=False, dest="specimen_metadata", required=False, help="specimen metadata tsv")
    parser.add_argument("-ex", "--experiment_metadata", default=False, dest="experiment_metadata", required=False, help="experiment metadata tsv")
    parser.add_argument("-rg", "--read_group_metadata", default=False, dest="read_group_metadata", required=False, help="read_group metadata tsv")
    parser.add_argument("-rm", "--relational_mapping", default=True, dest="relational_mapping", required=True, help="relational mapping json")
    parser.add_argument("-od", "--output_directory", default="output", dest="output_directory", required=False, help="output directory")
    args = parser.parse_args()

    main(args)