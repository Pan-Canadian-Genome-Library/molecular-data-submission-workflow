# Testing Guide

## Prerequisites

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

- **Test Environment**: These examples use the `cumulus_dev` profile for development testing in **OICR**. Please use `sd4h_dev` for any testing in **SD4H**.  
- **Authentication**: Replace `"test_token_here"` with your actual authentication token.
- **Network Requirements**: Ensure access to PCGL submission services in your testing environment.
- **Data Validation**: The provided test data is pre-validated; your actual data may require format adjustments.
- **Entity Registration**: Please make sure that `Study` and `Participant` entities are already registered for all test scenarios.
- **Data Model**: Please check the latest version of the [PCGL Base Data Model](https://drive.google.com/drive/u/1/folders/1vfNA7ajwh3WKkbVmswb6j9TuWKxaN9bB) to make sure your data conforms to the metadata requirements and dependencies.


## Running Tests

The repository includes comprehensive test datasets that you can use to verify the workflow functionality. Please be noted that the Submission Dependencies such as Study and Participant have already been registered beforehand. 
[TODO: add links to Study and Participant registration section when the full submission docs is available.]


### 🧪 **Test Scenario 1: No Pre-submitted Biospecimen Entities**
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

### 🧪 **Test Scenario 2: All Biospecimen Entities Pre-submitted**

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

### 🧪 **Test Scenario 3: Partial Biospecimen Entities Pre-submitted**

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

**What the mixed results indicate**:
- Some analyses succeeded with the provided metadata
- Others failed due to missing required metadata for their specific analysis type
- The exact success/failure split depends on your analysis types, metadata requirements, and biospecimen entity pre-submission state

## 📋 **Common Failure Issues**

When running the workflow, you may encounter these common issues across any test scenario:

### **Authentication and Network Issues**
- **Invalid authentication token**: Expired or incorrect authentication credentials
- **Network connectivity issues**: Problems connecting to PCGL submission services  
- **Service availability**: PCGL services temporarily unavailable

### **Biospecimen Entity Issues**
- **Biospecimen validation errors**: Issues with biospecimen metadata format or content
- **Entity dependency conflicts**: The provided biospecimen metadata conflict with existing entities
- **Missing metadata for required entities**: Some biospecimen entities expected but metadata not provided
- **Biospecimen entity not found**: Expected pre-submitted biospecimen entities don't exist in PCGL

### **Analysis and File Issues**  
- **Missing workflow metadata**: Required workflow metadata not provided for certain analysis types
- **Analysis type validation errors**: Analysis metadata doesn't meet requirements for the specified type
- **File metadata validation failures**: Issues with file metadata format or content
- **File transfer errors**: Issues uploading genomic files to object storage

### **System Issues**
- **Resource constraints**: Insufficient system resources during processing
- **Timeout errors**: Operations exceeding configured time limits

**👉 For detailed troubleshooting steps and solutions, see [Troubleshooting Guide](troubleshoot.md)**

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

**👉 See [Output Documentation](output.md) for comprehensive details about all output files and directories**


##  **Testing Best Practices**

1. **Start with Scenario 1 (No Pre-submitted Biospecimen)**: Verify full workflow functionality including biospecimen entity creation
2. **Progress to Scenario 2 (All Pre-submitted)**: Test file and analysis only submission pathway
3. **Test Scenario 3 (Partial Pre-submitted)**: Validate mixed entity state handling
4. **Review Batch Receipts**: Always check the receipt files after each test run to understand biospecimen entity processing
5. **Use Different Output Directories**: Keep test results organized by different scenario
6. **Check Work Directories**: Use work directory paths in receipts for detailed debugging information
7. **Understand Entity Dependencies**: Review which biospecimen entities your analyses require before testing


