# Testing Guide

## Prerequisites

Before testing the workflow, you can download the repository and its test data available locally.

## Getting the Repository

To get the full repository with all test data:

```bash
# Clone the repository
git clone https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow.git

# Navigate to the repository directory
cd molecular-data-submission-workflow

# Verify test data is available
ls tests/test_data/
```

## Running Tests

The repository includes comprehensive test datasets that you can use to verify the workflow functionality. Please be noted that the Submission Dependencies such as Study and Participant have already been registered beforehand. 


**Test with all provided metadata files (full test):**
```bash
nextflow run . \
    --study_id "TEST-CA" \
    --path_to_files_directory "tests/test_data/genomics" \
    --file_metadata "tests/test_data/analysis_meta/file_metadata.tsv" \
    --analysis_metadata "tests/test_data/analysis_meta/analysis_metadata.tsv" \
    --workflow_metadata "tests/test_data/analysis_meta/workflow_metadata.tsv" \
    --read_group_metadata "tests/test_data/biospecimen/read_group_metadata.tsv" \
    --experiment_metadata "tests/test_data/biospecimen/experiment_metadata.tsv" \
    --specimen_metadata "tests/test_data/biospecimen/specimen_metadata.tsv" \
    --sample_metadata "tests/test_data/biospecimen/sample_metadata.tsv" \
    --outdir test_results \
    -profile test,docker,cumulus_dev \
    --token "test_token_here"
```

**Test with minimal required files only:**
```bash
nextflow run . \
    --study_id "TEST-CA" \
    --path_to_files_directory "tests/test_data/genomics" \
    --file_metadata "tests/test_data/analysis_meta/file_metadata.tsv" \
    --analysis_metadata "tests/test_data/analysis_meta/analysis_metadata.tsv" \
    --outdir test_minimal \
    -profile test,docker,cumulus_dev \
    --token "test_token_here"
```

**Test with partial optional metadata:**
```bash
nextflow run . \
    --study_id "TEST-CA" \
    --path_to_files_directory "tests/test_data/genomics" \
    --file_metadata "tests/test_data/analysis_meta/file_metadata.tsv" \
    --analysis_metadata "tests/test_data/analysis_meta/analysis_metadata.tsv" \
    --experiment_metadata "tests/test_data/biospecimen/experiment_metadata.tsv" \
    --sample_metadata "tests/test_data/biospecimen/sample_metadata.tsv" \
    --outdir test_partial \
    -profile test,docker,cumulus_dev \
    --token "test_token_here"
    # workflow_metadata, read_group_metadata, specimen_metadata omitted
```
