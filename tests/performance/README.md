# Performance Testing Instructions

This document provides detailed instructions for conducting performance tests on the molecular data submission workflow using different data batch sizes and execution strategies.

## Test Data Overview

| Batch | File | Submissions | Data Size per Submission | Total Data Volume |
|-------|------|-------------|-------------------------|-------------------|
| 1 | `analysis_metadata_20G.tsv` | 20 | 20GB | 400GB |
| 2 | `analysis_metadata_39G.tsv` | 15 | 39GB | 585GB |
| 3 | `analysis_metadata_97G.tsv` | 10 | 97GB | 970GB |
| 4 | `analysis_metadata_194G.tsv` | 5 | 194GB | 970GB |
| 5 | `analysis_metadata_388G.tsv` | 3 | 388GB | 1.16TB |

## Prerequisites

1. **Environment Setup**
   ```bash
   # Ensure Nextflow is installed
   nextflow -version
   
   # Navigate to workflow directory
   cd /path/to/molecular-data-submission-workflow
   
   ```

2. **Data Preparation**
   - Ensure all batch metadata files are in `tests/performance/metadata/`
   - Verify corresponding data files are accessible at `../performance_test_data/`
   - Check sufficient disk space for trace files and outputs

## Test Scenario 1: Sequential Execution per Batch

Execute each batch independently with sequential submission processing.

### Commands

```bash
# Batch 1: 20GB submissions (sequential)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_20G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_20G \
    -profile test_sequential,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_sequential_20G.txt \
    -with-report tests/performance/trace_report/report_sequential_20G.html

# Batch 2: 39GB submissions (sequential)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_39G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_39G \
    -profile test_sequential,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_sequential_39G.txt \
    -with-report tests/performance/trace_report/report_sequential_39G.html

# Batch 3: 97GB submissions (sequential)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_97G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_97G \
    -profile test_sequential,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_sequential_97G.txt \
    -with-report tests/performance/trace_report/report_sequential_97G.html

# Batch 4: 194GB submissions (sequential)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_194G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_194G \
    -profile test_sequential,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_sequential_194G.txt \
    -with-report tests/performance/trace_report/report_sequential_194G.html

# Batch 5: 388GB submissions (sequential)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_388G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_388G \
    -profile test_sequential,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_sequential_388G.txt \
    -with-report tests/performance/trace_report/report_sequential_388G.html
```

## Test Scenario 2: Parallel Execution per Batch

Execute each batch independently allowing parallel submission processing.

### Commands

```bash
# Batch 1: 20GB submissions (parallel)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_20G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_parallel_20G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_parallel_20G.txt \
    -with-report tests/performance/trace_report/report_parallel_20G.html

# Batch 2: 39GB submissions (parallel)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_39G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_parallel_39G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_parallel_39G.txt \
    -with-report tests/performance/trace_report/report_parallel_39G.html

# Batch 3: 97GB submissions (parallel)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_97G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_parallel_97G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_parallel_97G.txt \
    -with-report tests/performance/trace_report/report_parallel_97G.html

# Batch 4: 194GB submissions (parallel)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_194G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_parallel_194G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_parallel_194G.txt \
    -with-report tests/performance/trace_report/report_parallel_194G.html

# Batch 5: 388GB submissions (parallel)
nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_388G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_parallel_388G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_parallel_388G.txt \
    -with-report tests/performance/trace_report/report_parallel_388G.html
```

## Test Scenario 3: Multi-VM Parallel Execution

Execute all 5 batches simultaneously from different VMs.

### VM1 Commands
```bash
# VM1: Execute Batch 1 & 2
nohup nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_20G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_vm1_20G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_vm1_20G.txt \
    -with-report tests/performance/trace_report/report_vm1_20G.html > vm1_20G.log 2>&1 &

nohup nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_39G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_vm1_39G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_vm1_39G.txt \
    -with-report tests/performance/trace_report/report_vm1_39G.html > vm1_39G.log 2>&1 &
```

### VM2 Commands
```bash
# VM2: Execute Batch 3 & 4
nohup nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_97G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_vm2_97G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_vm2_97G.txt \
    -with-report tests/performance/trace_report/report_vm2_97G.html > vm2_97G.log 2>&1 &

nohup nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_194G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_vm2_194G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_vm2_194G.txt \
    -with-report tests/performance/trace_report/report_vm2_194G.html > vm2_194G.log 2>&1 &
```

### VM3 Commands
```bash
# VM3: Execute Batch 5
nohup nextflow run . \
    --study_id "TEST-C3G" \
    --path_to_files_directory "../performance_test_data" \
    --file_metadata "tests/performance/metadata/file_metadata.tsv" \
    --analysis_metadata "tests/performance/metadata/analysis_metadata_388G.tsv" \
    --workflow_metadata "tests/performance/metadata/workflow_metadata.tsv" \
    --outdir test_performance_vm3_388G \
    -profile test_parallel,docker,sd4h_k8s_dev \
    --token "test_token_here" \
    -with-trace tests/performance/trace_report/trace_vm3_388G.txt \
    -with-report tests/performance/trace_report/report_vm3_388G.html > vm3_388G.log 2>&1 &
```

## Expected Outputs

After completing all tests, you should have both trace files (.txt) and execution report files (.html) in the same directory:

```
tests/performance/trace_report/
├── trace_sequential_20G.txt
├── report_sequential_20G.html
├── trace_sequential_39G.txt
├── report_sequential_39G.html
├── trace_sequential_97G.txt
├── report_sequential_97G.html
├── trace_sequential_194G.txt
├── report_sequential_194G.html
├── trace_sequential_388G.txt
├── report_sequential_388G.html
├── trace_parallel_20G.txt
├── report_parallel_20G.html
├── trace_parallel_39G.txt
├── report_parallel_39G.html
├── trace_parallel_97G.txt
├── report_parallel_97G.html
├── trace_parallel_194G.txt
├── report_parallel_194G.html
├── trace_parallel_388G.txt
├── report_parallel_388G.html
├── trace_vm1_20G.txt
├── report_vm1_20G.html
├── trace_vm1_39G.txt
├── report_vm1_39G.html
├── trace_vm2_97G.txt
├── report_vm2_97G.html
├── trace_vm2_194G.txt
├── report_vm2_194G.html
├── trace_vm3_388G.txt
└── report_vm3_388G.html
```

## Notes

- Ensure sufficient disk space before starting tests
- Monitor system resources to avoid overloading
- Sequential tests use `test_sequential` profile to limit parallelism
- Parallel tests use `test_parallel` profile for maximum throughput
- Replace `"test_token_here"` with actual authentication token
