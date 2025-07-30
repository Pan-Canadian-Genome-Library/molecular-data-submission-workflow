#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECK_CLINICAL } from '../../../modules/local/checkclinical/'

workflow test_check_input {

    //params.read_group_metadata == null ? [[]] : [params.read_group_metadata]
    CHECK_CLINICAL(
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
}
