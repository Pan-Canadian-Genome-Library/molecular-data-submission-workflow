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
[TODO: add links to Study and Participant registration section when the full submission docs is available.]


### 🧪 **Test Scenario 1: No Pre-submitted Biospecimen Entities**
**Purpose**: Test the workflow when we assume no biospecimen entities (samples, specimens, experiments, read groups) have been submitted to PCGL in advance. This scenario requires the workflow to create and submit all biospecimen entities during submission.

**What this tests**: 
- Complete biospecimen entity creation
- Full metadata validation and submission of biospecimen data
- Clinical data submission process for all biospecimen data
- End-to-end submission process including biospecimen submission

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

**Expected Terminal Output**: Upon successful completion, you will see a comprehensive summary box:

```
╔══════════════════════════════════════════════════════════════════════════════════╗
║                       🎉 WORKFLOW COMPLETED! 🎉
╠══════════════════════════════════════════════════════════════════════════════════╣
║  Study: TEST-CA
║  Batch ID: batch_20250829_130026
║  Total in this batch:    5
║  ✅ Successful submissions: 5  
║  ❌ Failed submissions:     0
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

### 🧪 **Test Scenario 2: All Biospecimen Entities Pre-submitted**

**Purpose**: Test the workflow when we assume all biospecimen entities (samples, specimens, experiments, read groups) have already been submitted to PCGL in advance. This scenario only requires file and analysis metadata for submission.

**What this tests**:
- Analysis-only submission pathway (biospecimen entities already exist)
- Minimal metadata validation (no biospecimen data required)
- Workflow behavior when biospecimen dependencies are satisfied or not
- Error handling for missing optional metadata

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

**Expected Terminal Output**: Upon successful completion, you will see a comprehensive summary box:

```
╔══════════════════════════════════════════════════════════════════════════════════╗
║                       🎉 WORKFLOW COMPLETED! 🎉
╠══════════════════════════════════════════════════════════════════════════════════╣
║  Study: TEST-CA
║  Batch ID: batch_20250829_130026
║  Total in this batch:    5
║  ✅ Successful submissions: 0
║  ❌ Failed submissions:     5
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

**Note**: This test may result in failures if your analysis type requires workflow metadata. Check the batch receipt for specific error details such as:
- `"Error: workflow-meta is required for analysis type 'sequenceAlignment'"`
- `"Missing required biospecimen metadata dependencies"`

### 🧪 **Test Scenario 3: Partial Biospecimen Entities Pre-submitted**

**Purpose**: Test the workflow when some biospecimen entities have been submitted to PCGL in advance, but others have not. This tests the workflow's ability to handle mixed biospecimen entity states.

**What this tests**:
- Mixed biospecimen entity submission scenarios
- Selective biospecimen metadata processing
- Conditional validation logic for existing vs. new entities
- Workflow robustness with partial biospecimen dependencies
- Hybrid submission pathway (some entities exist, others need creation)
- Robustness with incomplete optional data

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

**Expected Terminal Output**: Upon successful completion, you will see a comprehensive summary box:

```
╔══════════════════════════════════════════════════════════════════════════════════╗
║                       🎉 WORKFLOW COMPLETED! 🎉
╠══════════════════════════════════════════════════════════════════════════════════╣
║  Study: TEST-CA
║  Batch ID: batch_20250829_140132
║  Total in this batch:    5
║  ✅ Successful submissions: 2
║  ❌ Failed submissions:     3
║
║  📋 BATCH RECEIPT GENERATED:
║  📁 JSON RECEIPT Location: /path/to/work/dir/20250829_140132_batch_receipt.json
║  📁 TSV RECEIPT Location: /path/to/work/dir/20250829_140132_batch_receipt.tsv
║
║  ℹ️  The batch receipt contains:
║     • Summary of all processed analyses
║     • Status of each submission step
║     • Upload results and analysis IDs
║     • Error details for any failed analyses
╚══════════════════════════════════════════════════════════════════════════════════╝
```

**What the mixed results indicate**:
- Some analyses succeeded with the provided partial metadata
- Others failed due to missing required metadata for their specific analysis type
- The exact success/failure split depends on your analysis types and metadata requirements

## 📋 **Understanding Test Results Output**

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


## 🔍 **Common Test Scenarios and Expected Outcomes**

| Test Scenario | Expected Outcome | Common Issues |
|---------------|------------------|---------------|
| **No Pre-submitted Biospecimen** | Mixed SUCCESS/FAILED | Network issues, biospecimen validation errors, invalid token issues|
| **All Pre-submitted Biospecimen** | Mixed SUCCESS/FAILED | Missing workflow_metadata, invalid token issues |
| **Partial Pre-submitted Biospecimen** | Mixed SUCCESS/FAILED | Entity dependency conflicts, missing metadata issues, invalid token issues |

## 💡 **Testing Best Practices**

1. **Start with Scenario 1 (No Pre-submitted Biospecimen)**: Verify full workflow functionality including biospecimen entity creation
2. **Progress to Scenario 2 (All Pre-submitted)**: Test file and analysis only submission pathway
3. **Test Scenario 3 (Partial Pre-submitted)**: Validate mixed entity state handling
4. **Review Batch Receipts**: Always check the receipt files after each test run to understand biospecimen entity processing
5. **Use Different Output Directories**: Keep test results organized by different scenario
6. **Check Work Directories**: Use work directory paths in receipts for detailed debugging information
7. **Understand Entity Dependencies**: Review which biospecimen entities your analyses require before testing

## ⚠️ **Important Notes**

- **Test Environment**: These examples use the `cumulus_dev` profile for development testing in OICR. Please use `sd4h_dev` for any testing in SD4H.  
- **Authentication**: Replace `"test_token_here"` with your actual authentication token
- **Network Requirements**: Ensure access to PCGL submission services
- **Data Validation**: Test data is pre-validated; your actual data may require format adjustments
- **Entity Registration**: Study and Participant entities are already registered for test scenarios
