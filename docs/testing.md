## Testing the Workflow

The repository includes comprehensive test datasets that you can use to verify the workflow functionality:

**Test with all provided metadata files (full test):**
```bash
nextflow run . \
    --study_id "TEST-CA" \
    --token "test_token_here" \
    --path_to_files_directory "tests/test_data/genomics" \
    --file_metadata "tests/test_data/analysis_meta/file_metadata.tsv" \
    --analysis_metadata "tests/test_data/analysis_meta/analysis_metadata.tsv" \
    --workflow_metadata "tests/test_data/analysis_meta/workflow_metadata.tsv" \
    --read_group_metadata "tests/test_data/biospecimen/read_group_metadata.tsv" \
    --experiment_metadata "tests/test_data/biospecimen/experiment_metadata.tsv" \
    --specimen_metadata "tests/test_data/biospecimen/specimen_metadata.tsv" \
    --sample_metadata "tests/test_data/biospecimen/sample_metadata.tsv" \
    --skip_upload true \
    --outdir test_results \
    -profile test,docker,cumulus_dev
```

**Test with minimal required files only:**
```bash
nextflow run . \
    --study_id "TEST-CA" \
    --token "test_token_here" \
    --path_to_files_directory "tests/test_data/genomics" \
    --file_metadata "tests/test_data/analysis_meta/file_metadata.tsv" \
    --analysis_metadata "tests/test_data/analysis_meta/analysis_metadata.tsv" \
    --skip_upload true \
    --outdir test_minimal \
    -profile test,docker,cumulus_dev
```

**Test with partial optional metadata:**
```bash
nextflow run . \
    --study_id "TEST-CA" \
    --token "test_token_here" \
    --path_to_files_directory "tests/test_data/genomics" \
    --file_metadata "tests/test_data/analysis_meta/file_metadata.tsv" \
    --analysis_metadata "tests/test_data/analysis_meta/analysis_metadata.tsv" \
    --experiment_metadata "tests/test_data/biospecimen/experiment_metadata.tsv" \
    --sample_metadata "tests/test_data/biospecimen/sample_metadata.tsv" \
    --skip_upload true \
    --outdir test_partial \
    -profile test,docker,cumulus_dev
    # workflow_metadata, read_group_metadata, specimen_metadata omitted
```
