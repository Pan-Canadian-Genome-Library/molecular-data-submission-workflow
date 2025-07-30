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
        params.data_directory == null ? [[]] : [params.data_directory]
    )

    CHECK_SUBMISSION_DEPENDENCIES.out.analysis_channels.subscribe{
        println "${it.meta.id} META:${it.meta}\n${it.meta.id} ANALYSIS:${it.analysis}\n${it.meta.id} CLINICAL:${it.clinical}\n${it.meta.id} FILES:${it.files}\n${it.meta.id} STATUS:${it.status}\n${it.meta.id} STATUS_FILE:${it.status_file}"}


}
