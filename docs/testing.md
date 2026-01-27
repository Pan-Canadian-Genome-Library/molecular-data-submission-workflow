# PCGL Molecular Data Submission Workflow: Testing Guide

## 📋 **Prerequisites**

Before testing the workflow, ensure you have the proper environment set up and can download the repository with its test data available locally.

### Environment Setup
For detailed instructions on installing Nextflow, container engines, and obtaining API tokens, see [Environment Setup](../README.md#environment-setup) in the main README.

### Getting the Repository

To get the full repository with all test data:

```bash
# Clone the repository
git clone https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow.git

# Navigate to the repository directory
cd molecular-data-submission-workflow

# Verify test data is available
ls tests/test_data/
```

## ⚠️ **Important Notes**

- **Test Environment**: 
  - These examples use the `cumulus_dev` profile for testing in **OICR Dev**. 
  - Please use `sd4h_dev` for any testing in **SD4H Dev**, 
  - Please use `sd4h_staging` for any testing in **SD4H Staging**.  
- **Authentication**: Replace `"test_token_here"` with your actual authentication token.
- **Network Requirements**: Ensure access to PCGL submission services in your testing environment.
- **Data Model**: Please check the latest version of the [PCGL Base Data Model](https://drive.google.com/drive/u/1/folders/1vfNA7ajwh3WKkbVmswb6j9TuWKxaN9bB) to ensure your data conforms to the metadata requirements and dependencies.
- **Entity Registration**: Please make sure that `Study` and `Participant` entities are already registered for all test scenarios.
- **Test Dataset**: 
  - The [provided test data](../tests/test_data/) is pre-formatted and compliant with the latest PCGL Base Data Model. Please remember to replace **studyId** in `tests/test_data/analysis_meta/analysis_metadata.tsv` to **your_study_id**. 
  - When using your own data, refer to the [Input Documentation](input.md) for formatting requirements and data preparation guidelines.


## 🧪 **Running Tests**

The repository includes comprehensive test datasets that you can use to verify the workflow functionality. Please be noted that the Submission Dependencies such as Study and Participant have already been registered beforehand. 
[**TODO**: add links to Study and Participant registration section when the full submission docs is available.]


### **Test Scenario 1: No Pre-submitted Biospecimen Entities**
**Purpose**: Test the workflow when we assume no biospecimen entities (samples, specimens, experiments, read groups) have been submitted to PCGL in advance. This scenario requires the workflow to create and submit all biospecimen entities during submission.

**What this tests**: 
- Complete biospecimen entity creation and submission
- Full metadata validation for all biospecimen data
- End-to-end submission process including biospecimen submission to clinical system, analysis and file metadata submission to file-manager, and genomic files submission to object storage

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

### **Test Scenario 2: All Biospecimen Entities Pre-submitted**

**Purpose**: Test the workflow when we assume all biospecimen entities (samples, specimens, experiments, read groups) have already been submitted to PCGL in advance. This scenario only requires file and analysis metadata for submission.

**What this tests**:
- Workflow's behavior when biospecimen dependencies are already satisfied
- Validation against existing biospecimen metadata
- Error handling for missing optional metadata
- Analysis-only submission pathway when biospecimen entities already exist in PCGL clinical system, analysis and file metadata submission to file-manager, and genomic files submission to object storage

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

### **Test Scenario 3: Partial Biospecimen Entities Pre-submitted**

**Purpose**: Test the workflow when we assume some biospecimen entities have been submitted to PCGL in advance, but others have not. This tests the workflow's ability to handle mixed biospecimen entity states.

**What this tests**:
- Mixed biospecimen entity states (some pre-submitted, some new)
- Workflow's ability to handle partial metadata scenarios and entity dependencies
- Selective biospecimen metadata processing based on existing entity states
- Hybrid submission pathway combining existing and new entity submission to clinical system, analysis and file metadata submission to file-manager, and genomic files submission to object storage

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

## 📊 **Expected Terminal Output**

When you run any test scenario, the workflow will display input parameters, execute pipeline processes in real-time, and conclude with either a successful completion summary or an early termination message if critical issues are detected.

### **Input Parameters Display**

At the start of workflow execution, you'll see a summary of all input parameters, for example:

```
🔧 Input parameters:
   - study_id: TEST-CA
   - analysis_metadata: tests/test_data/analysis_meta/analysis_metadata.tsv
   - file_metadata: tests/test_data/analysis_meta/file_metadata.tsv
   - workflow_metadata: tests/test_data/analysis_meta/workflow_metadata.tsv
   - read_group_metadata: tests/test_data/biospecimen/read_group_metadata.tsv
   - experiment_metadata: tests/test_data/biospecimen/experiment_metadata.tsv
   - specimen_metadata: tests/test_data/biospecimen/specimen_metadata.tsv
   - sample_metadata: tests/test_data/biospecimen/sample_metadata.tsv
   - path_to_files_directory: tests/test_data/genomics
   - skip_upload: false
   - allow_duplicates: true
```

**What this shows**:
- **Study identification**: Confirms which study is being processed
- **Configuration confirmation**: Verifies the workflow received all your input parameters correctly
- **File paths validation**: Shows the exact paths to metadata files being used
- **Parameter settings**: Displays key workflow behavior settings like `skip_upload` and `allow_duplicates`


### **Pipeline Process Execution**

After the input parameters, you'll see Nextflow executing individual pipeline processes in real-time, for example:

```
[0a/132c7b] PCG…CK_DEPENDENCIES (TEST-CA) | 1 of 1 ✔
[32/3ff298] PCG…:ANALYSIS_SPLIT (TEST-CA) | 1 of 1 ✔
[e2/669398] PCG…E_CLINICAL (analysis_005) | 5 of 5 ✔
[c9/d27988] PCG…SUBMISSION (analysis_005) | 5 of 5 ✔
[1b/a88c56] PCG…D_GENERATE (analysis_005) | 4 of 4 ✔
[62/43ef84] PCG…D_VALIDATE (analysis_005) | 4 of 4 ✔
[d4/0eda2e] PCG…N_METADATA (analysis_005) | 4 of 4 ✔
......
```

**What this shows**:
- **Process IDs**: Each line starts with a unique process execution ID (e.g., `[0a/132c7b]`)
- **Process names**: Abbreviated process names (e.g., `PCG…CK_DEPENDENCIES` for dependency checking)
- **Input context**: What the process is operating on (e.g., study ID, analysis ID, etc.)
- **Execution progress**: Number of tasks completed vs. total tasks (e.g., `5 of 5`)
- **Status indicators**: ✔ indicates successful completion of that process
- **Real-time updates**: Processes appear as they execute, showing workflow progress

### **Case 1: Workflow Stops Due to Minimum Requirements Not Met**

If the workflow detects that minimum requirements for data submission are not met, it will stop early and display:

```
╔══════════════════════════════════════════════════════════════════════════════╗
║                        🚨 WORKFLOW STOPPED                                   ║
║              Minimum requirements for data submission not met!               ║
╚══════════════════════════════════════════════════════════════════════════════╝

❌ Study: TEST-CA

🔍 Issues found:
--------------------------------------------------------------------------------
• Missing required metadata file: workflow_metadata.tsv
• Invalid file path in file_metadata.tsv: sample001.bam (file not found)
• Study ID mismatch: Expected 'TEST-CA' but found 'TEST-CB' in analysis_metadata.tsv
--------------------------------------------------------------------------------

📋 Common solutions:
   • Ensure that access token with correct submission scope is provided
   • Ensure that the study has been registered
   • Ensure all required metadata files are provided and correctly formatted
   • Check that file paths in metadata files are valid and files exist
   • Verify that all required columns are present in metadata files
   • Ensure study_id matches across all metadata files

💡 For detailed error information, check: /path/to/work/dir/status.yml

Please fix the above issues and re-run the workflow.
```

**What this output means**:
- **🚨 WORKFLOW STOPPED**: The workflow has terminated early before processing any analyses
- **Study**: Shows which study the workflow was attempting to process
- **Issues found**: Detailed list of specific problems that prevented the workflow from starting
- **Common solutions**: Quick reference guide for resolving typical configuration issues
- **Status file location**: Path to detailed error information for further investigation

❌ **When you see this output**: The workflow has detected critical issues that prevent data submission from proceeding. No analyses will be processed until these fundamental requirements are resolved. Fix the listed issues and re-run the workflow.

### **Case 2: Workflow Completes Successfully**

Upon successful completion of any test scenario, you will see a comprehensive summary box with mixed SUCCESS/FAILED results:

```
╔══════════════════════════════════════════════════════════════════════════════════╗
║                       🎉 WORKFLOW COMPLETED! 🎉
╠══════════════════════════════════════════════════════════════════════════════════╣
║  Study: TEST-CA
║  Batch ID: batch_20250829_130026
║  Total in this batch:    5
║  ✅ Successful submissions: 2  
║  ❌ Failed submissions:     3
║
║  📋 BATCH RECEIPT GENERATED:
║  📁 JSON RECEIPT Location: /path/to/work/dir/20250829_130026_batch_receipt.json
║  📁 TSV RECEIPT Location: /path/to/work/dir/20250829_130026_batch_receipt.tsv
║
║  ℹ️  The batch receipt contains:
║     • Summary of all processed analyses
║     • Status of each submission step
║     • Upload results and analysis IDs
║     • Error details for any failed analyses
╚══════════════════════════════════════════════════════════════════════════════════╝
```

**What this output means**:
- **Study**: Confirms the study being processed
- **Batch ID**: Unique identifier for this submission batch  
- **Success/Failure counts**: Quick overview of submission results
- **Receipt locations**: Direct paths to detailed batch receipt files
- **Summary information**: What you'll find in the receipt files

**What the results indicate**:
- **Successful analyses**: These analyses had all required metadata and passed validation successfully.
- **Failed analyses**: These encountered issues during processing and require attention before resubmission.
- **Result variability**: The success/failure depends on your analysis types, metadata completeness, and biospecimen entity pre-submission state.
- **Next steps**: Review the batch receipt to identify specific failure reasons and correct the issues before re-running the workflow.


## 📋 **Understanding Test Results**

### **Finding Your Batch Receipt**

The workflow completion summary shows the exact location of your batch receipt files. You can find them in two ways:

**Method 1: From Terminal Output**
The completion summary provides direct paths to your receipt files:
```
║  📁 JSON RECEIPT Location: /path/to/work/dir/20250829_130026_batch_receipt.json
║  📁 TSV RECEIPT Location: /path/to/work/dir/20250829_130026_batch_receipt.tsv
```

Simply copy these paths to access your receipt files directly.

**Method 2: Standard Output Directory Structure**
Receipts are also copied to your specified output directory:
```
<outdir>/
├── receipt_aggregate/
│   ├── YYYYMMDD_HHMMSS_batch_receipt.json    # Detailed JSON format
│   ├── YYYYMMDD_HHMMSS_batch_receipt.tsv     # Tabular format
```

### **Understanding Your Receipt**

For comprehensive information about interpreting your batch receipt, including:
- Detailed field descriptions
- Error message meanings
- Troubleshooting guidance
- Best practices for receipt analysis

**👉 See [Receipt Documentation](receipt.md) for complete details**

### **Examining the Specified Output Directory**

Beyond the batch receipt files, your specified output directory (`--outdir`) contains additional important information about the workflow execution:

**Complete Output Structure**:
Your output directory will contain various subdirectories with intermediate files, logs, and results from each workflow process. This includes:
- Process-specific input/output files
- Intermediate data files
- Process execution logs
- Pipeline execution reports

**👉 See [Output Documentation](output.md) [**TODO**] for comprehensive details about all output files and directories**

## 📋 **Common Failure Issues**

When running the workflow, you may encounter these common issues across any test scenario:

### **Authentication and Network Issues**
- **Invalid authentication token**: Expired or incorrect authentication credentials will result in authentication failures (typically HTTP `403` Forbidden errors).
- **Network connectivity issues**: Problems connecting to PCGL submission services  
- **Service availability**: PCGL services temporarily unavailable

### **Submission Dependencies Issues**
- **Unregistered study**: Study ID not found in PCGL system - contact PCGL Administrator to register your study before submission
- **Unregistered participants**: One or more participant IDs in your submission batch are not registered in the PCGL clinical system
- **Biospecimen metadata validation errors**: Biospecimen metadata files contain format errors, missing required fields, or invalid data values
- **Biospecimen entity dependency conflicts**: Provided biospecimen metadata conflicts with existing entities already registered in PCGL
- **Missing biospecimen metadata**: Required biospecimen entities (samples, specimens, experiments, read groups) are expected but metadata files not provided
- **Missing biospecimen entities**: Expected pre-submitted biospecimen entities do not exist in PCGL clinical system - verify entity IDs and registration status

### **Analysis and File Issues**  
- **Missing workflow metadata**: Required `workflow_metadata.tsv` file not provided for analysis types that mandate workflow provenance information (e.g., variant calling pipelines, assembly workflows)
- **Analysis type validation errors**: Analysis metadata values fail schema validation for the specified analysis type - check required fields, controlled vocabulary terms, and data type constraints in `analysis_metadata.tsv`
- **File validation failures**: Genomic file metadata contains formatting errors, missing mandatory fields, or invalid values; alternatively, genomic files referenced in `file_metadata.tsv` are missing from the specified directory path or fail checksum verification
- **File transfer errors**: Genomic file upload to PCGL object storage fails due to network issues, authentication problems, or object storage service unavailability

### **System Issues**
- **Resource constraints**: Insufficient system resources during processing
- **Timeout errors**: Operations exceeding configured time limits

**👉 For detailed troubleshooting steps and solutions, see [Troubleshooting Guide](troubleshooting.md)** [**TODO**]



## 💡 **Testing Best Practices**
- **Understand Submission Dependencies**: Review entity hierarchy (Study → Participant → Sample → Specimen → Experiment → Read Group) and analysis requirements before test execution
- **Start with Scenario 1 (No Pre-submitted Biospecimen)**: Execute complete end-to-end workflow testing including biospecimen entity creation, validation, and submission to verify full pipeline functionality
- **Progress to Scenario 2 (All Pre-submitted)**: Test analysis-only submission pathway when biospecimen dependencies are satisfied - validates workflow behavior with existing entity references
- **Test Scenario 3 (Partial Pre-submitted)**: Validate hybrid submission mode handling mixed biospecimen entity states and selective metadata processing
- **Nextflow Process Status Interpretation**: Nextflow processes show completion status due to error handling mechanisms. Therefore failed processes still appear "completed" but with error exit codes captured in batch receipt
- **Review Batch Receipts**: Examine JSON or TSV receipt files after each test execution to analyze submission status and failure diagnostics. Failed processes trigger automatic skipping of dependent downstream processes. Please analyze the process entries which were captured in chronological execution order in the JSON receipt to identify the root cause of the failed submission.
- **Use Different Output Directories**: Maintain separate `--outdir` paths for each test scenario to prevent result contamination and enable parallel testing
- **Check Work Directories**: Utilize work directory paths documented in batch receipts for process-level debugging and intermediate file inspection




