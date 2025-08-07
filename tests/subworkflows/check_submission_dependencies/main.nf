#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECK_SUBMISSION_DEPENDENCIES } from '../../../subworkflows/local/check_submission_dependencies.nf'

workflow test_check_input {

    //params.read_group_metadata == null ? [[]] : [params.read_group_metadata]
    CHECK_SUBMISSION_DEPENDENCIES(
        [params.study_id],
        [params.token],
        [params.file_metadata],
        [params.analysis_metadata],
        params.workflow_metadata == null ? [[]] : [params.workflow_metadata],
        params.read_group_metadata == null ? [[]] : [params.read_group_metadata],
        params.experiment_metadata == null ? [[]] : [params.experiment_metadata],
        params.specimen_metadata == null ? [[]] : [params.specimen_metadata],
        params.sample_metadata == null ? [[]] : [params.sample_metadata],
        params.data_directory
    )

    CHECK_SUBMISSION_DEPENDENCIES.out.analysis_channels.subscribe{
    println """
    ${it.meta.id}\tMETA
    \tID : ${it.meta.id}
    \tstudy : ${it.meta.study}
    \ttype : ${it.meta.type}
    ${it.meta.id}\tANALYSIS
    \tanalysis : ${it.analysis.analysis}
    \tworkflow : ${it.analysis.workflow}
    \tfiles : ${it.analysis.files}
    ${it.meta.id}\tCLINICAL
    \tspecimen : ${it.clinical.specimen}
    \tsample : ${it.clinical.sample}
    \texperiment : ${it.clinical.experiment}
    \tread_group : ${it.clinical.read_group}
    ${it.meta.id}\tFILES
    \tfiles : \n${"\t\t"+it.files.join("\n\t\t")}
    ${it.meta.id}\tSTATUS_FILE
    \tstatus_file : ${it.status_file}
    """
        }


}
