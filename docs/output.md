# Pan-Canadian-Genome-Library/molecular-data-submission-workflow: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the output directory after the pipeline has finished. All paths are relative to the top-level results directory.

The results directory can be specified via `--outdir` otherwise the default will be `out`.

For tracking, the suggestion folder structure for output is :
```
out
|--batch1
|--batch2
```
such that the supplied output directory is `--outdir out/batch1` and `--outdir out/batch2` per instance of workflow. This ensures receipts are written to folders specific to the batch.

## Pipeline overview

The output contents of the pipeline is organized per process and organized according by chronological order in pipeline: 

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:
- [check_dependencies](#check_dependencies) - preflight check validating servers, studies, and input TSVs
- [analysis_split](#analysis_split) - runs workflow per record of input `analysis.tsv`
- [validate_clinical](#validate_clinical) - checks biospecimen content for discrepancies and dependencies
- [clinical_submission](#clinical_submission) - submits biospecimen data to clinical-submission service
- [payload_generate](#payload_generate) - aggregates analysis, file and workflow metadata into JSON payload
- [payload_validate](#payload_validate) - validates payload against file manager schema
- [validation_metadata](#validation_metadata) - validates file data with read group and payload
- [validation_crosscheck](#validation_crosscheck) - check md5 between payload and files
- [seqkit_seq](#seqkit_seq) - validates compressed/uncompressed FASTQ per file
- [samtools_quickcheck](#samtools_quickcheck) - validates BAM/CRAM format per file
- [bcftools_view](#bcftools_view) - validates VCF/BCF format per file
- [song_submit](#song_submit) - submits JSON payload to File manager to reserve `objectIDs`
- [song_manifest](#song_manifest) - generates a TXT to identify local file paths
- [score_upload](#score_upload) - upload data on File Transfer
- [song_publish](#song_publish) - uploaded data is made live and discoverable
- [song_getanalysis](#song_getanalysis) - return generated analysisId
- [receipt_generate](#receipt_generate) - combines `status.yml` from previous steps into a single `status.yml` per analyis
- [receipt_aggregate](#receipt_aggregate) - combines `status.yml` from all analyses into a single batch `status.yml`
- [pipeline_info](#pipeline_info) - Reports metrics 

## Workflow steps
### check_dependencies
- `{studyId}_pcgl_molecular_data_submission_workflow_check_submission_dependencies_check_dependencies_status.yml`
   - occurs per `Study`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### analysis_split
- `{studyId}/{submitterAnalysisId}/{submitterAnalysisId}_check_submission_dependencies_analysis_split_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### validate_clinical
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_check_submission_dependencies_validate_clinical_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### clinical_submission
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_check_submission_dependencies_clinical_submission_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### payload_generate
- `{submitterAnalysisId}_payload.json`
   - occurs per `Analysis`
   - for example see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/payload/payload_bam.json)
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_metadata_payload_generation_payload_generate_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### payload_validate
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_metadata_payload_generation_payload_validate_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### validation_metadata
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_validation_validation_metadata_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### validation_crosscheck
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_validation_validation_crosscheck_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### seqkit_seq
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_validation_file_integrity_seqkit_seq_{fileName}_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### bcftools_view
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_validation_file_integrity_bcftools_view_{fileName}_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### samtools_quickcheck
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_validation_file_integrity_samtools_quickcheck_{fileName}_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### song_submit
- `{submitterAnalysisId}_analysis_id.txt`
   - occurs per `Analysis`
   - a text file containing the corresponding SONG `analysisId`
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_upload_song_submit_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### song_manifest
- `{submitterAnalysisId}_manifest.txt`
   - occurs per `Analysis`
   - a table containing `analysisId` header followed by a body of : `objectId`,`fileName` and `md5Sum`
   - example:

```
792b1f90-8775-4e9b-ab1f-9087755e9bd6
f27e60ae-b25e-537e-947b-7bed6246a963    C0HVY.2_r1.fq.gz       64cf635dbc54f53cae2cb03ec8e8471b
2ad85792-7a8c-5d50-bf26-e567c36f8e9d    C0HVY.2_r2.fq.gz       ceb7e66d031cb894fa9d1d3f8da65fc7
```

- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_upload_song_manifest_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### score_upload
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_upload_score_upload_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### song_publish
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_upload_song_publish_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### song_getanalysis
- `{submitterAnalysisId}_pcgl_molecular_data_submission_workflow_data_upload_song_getanalysis_status.yml`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
- `{submitterAnalysisId}_{analysisId}.analysis.json`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#3-process-level)
   - see [example](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/tests/test_data/status/test_status1.yml)
### receipt_generate
- `{submitterAnalysisId}_receipt.json`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#2-analysis-level)
### receipt_aggregate
- `{timestamp}_batch_receipt.json`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#2-analysis-level) breakdown and [JSON format](https:https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#json-format-_batch_receiptjson)
- `{timestamp}_batch_receipt.tsv`
   - occurs per `Analysis`
   - for contents breakdown see : [URL](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow/blob/main/docs/receipt.md#json-format-_batch_receiptjson)
### pipeline_info
- `pipeline_software_versions.yml`
  - `yml` file containing the software version of each process
- `pipeline_dag_{timestamp}.html`
  - directed acylic graph of nextflow processes, input and output
- `execution_trace_{timestamp}.txt`
  - TSV summary of processes wih nextflow task IDs, duration, resource usage, and summary status
- `execution_timeline_{timestamp}.html`
  - timeline chart of each the pipeline breaking down duration of each process and I/O speed
- `execution_report_{timestamp}.html`
  - html file summary of the workflow with the command run, resource allocation, job duration and task status summary



