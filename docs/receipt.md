# PCGL Molecular Data Submission Workflow: Understanding Batch Receipt

The batch receipt is a comprehensive execution report generated after processing all analyses in your submission batch. It provides detailed status information about the success or failure of each analysis across three hierarchical levels:

- **[Batch Level](#1-batch-level)**: Aggregated submission statistics for the entire batch which contains multiple analyses and provides overall submission summary and tracking information during batch processing
- **[Analysis Level](#2-analysis-level)**: Individual analysis submission results which contains analysis-specific metadata, PCGL system identifiers, and published status  
- **[Process Level](#3-process-level)**: Granular workflow process execution details which tracks each Nextflow task/step within the pipeline with process name, execution status, exit codes, timestamps, and debugging information

## 📋 **Receipt File Formats**

The workflow generates receipt files in two formats:

### **JSON Format (`*_batch_receipt.json`)**
- **Purpose**: Contains complete detailed information in a structured JSON format
- **Use Case**: Programmatic processing, detailed troubleshooting, comprehensive audit trails
- **Content**: Full hierarchical data including all process details and error messages

### **TSV Format (`*_batch_receipt.tsv`)**
- **Purpose**: Tabular format for easy viewing and analysis
- **Use Case**: Quick overview, summary reporting
- **Content**: Flattened summary data

## 📊 **JSON Format Receipt Structure Hierarchy**

The receipt information is organized in three levels:

### **1. Batch Level** 
Contains overall batch submission statistics:

| Field | Description | Example |
|-------|-------------|---------|
| `batch_id` | Unique identifier for this batch submission | `"20250828_165136"` |
| `generated_at` | Timestamp when the receipt was generated | `"2025-08-28T20:51:38.068683Z"` |
| `total_analyses` | Total number of analyses processed | `5` |
| `successful_analyses` | Number of analyses that completed successfully | `2` |
| `failed_analyses` | Number of analyses that failed | `3` |

### **2. Analysis Level**
For each analysis in the batch:

| Field | Description | Example |
|-------|-------------|---------|
| `submitter_analysis_id` | Your unique identifier for the analysis | `"analysis_004"` |
| `file_manager_analysis_id` | PCGL system assigned analysis ID (if successful) | `"0fe74e99-f30e-40b3-a74e-99f30eb0b302"` or `"Not applicable"` |
| `overall_status` | Final status of the entire analysis | `"SUCCESS"` or `"FAILED"` |
| `analysis_type` | Type of genomic analysis performed | `"sequenceAlignment"`, `"variantCall"` |
| `study_id` | Study identifier | `"TEST-CA"` |
| `analysis_state` | Current state in PCGL system (if successful)| `"PUBLISHED"` or `"Not applicable"` |
| `published_at` | Published timestamp (if successful) | `"2025-08-28T20:52:00Z"` or `"Not applicable"` |
| `generated_at` | When this analysis receipt was generated | `"2025-08-28T20:51:30.966475Z"` |

### **3. Process Level**
Each analysis contains detailed information about individual workflow processes:

| Field | Description | Example |
|-------|-------------|---------|
| `process` | Full process name from the workflow | `"PCGL:MOLECULAR_DATA_SUBMISSION_WORKFLOW:METADATA_PAYLOAD_GENERATION:PAYLOAD_GENERATE"` |
| `status` | Process execution status | `"SUCCESS"`, `"FAILED"` |
| `exit_code` | System exit code (0=success, 1=failed) | `0`, `1` |
| `timestamp` | When the process completed | `"2025-08-28T20:51:18+00:00"` |
| `work_directory` | Execution directory path for debugging | `"/path/to/work/dir"` |
| `details` | Process-specific information and parameters | `{"analysis_id": "analysis_004", "error_message": "",...}` |
| `details/error_message` | Detailed error description (if failed) | `"Error: workflow-meta is required for analysis type 'sequenceAlignment'"` |

## 📝 **JSON Receipt Example**

```json
{
  "batch_id": "20250828_165136",
  "generated_at": "2025-08-28T20:51:38.068683Z",
  "total_analyses": 5,
  "successful_analyses": 0,
  "failed_analyses": 5,
  "analyses": [
    {
      "submitter_analysis_id": "analysis_004",
      "file_manager_analysis_id": "Not applicable",
      "overall_status": "FAILED",
      "analysis_type": "sequenceAlignment",
      "study_id": "TEST-CA",
      "analysis_state": "Not applicable",
      "published_at": "Not applicable",
      "generated_at": "2025-08-28T20:51:30.966475Z",
      "processes": [
        {
          "process": "CHECK_SUBMISSION_DEPENDENCIES:ANALYSIS_SPLIT",
          "status": "SUCCESS",
          "exit_code": 0,
          "timestamp": "2025-08-28 20:51:01.287276",
          "work_directory": "/work/8b/39a4c2003dc7eacffe00b1bfb8587e",
          "details": {
            "analysis_id": "analysis_004",
            "error_message": ""
          }
        },
        {
          "process": "METADATA_PAYLOAD_GENERATION:PAYLOAD_GENERATE",
          "status": "FAILED",
          "exit_code": 1,
          "timestamp": "2025-08-28T20:51:18+00:00",
          "work_directory": "/work/a6/7713b639c9fb893ab0e0dc684ac263",
          "details": {
            "analysis_id": "analysis_004",
            "payload_file": "analysis_004_payload.json",
            "file_meta": "files.tsv",
            "analysis_meta": "analysis.tsv",
            "workflow_meta": "",
            "exit_on_error_enabled": "false",
            "error_message": "Error: workflow-meta is required for analysis type 'sequenceAlignment'"
          }
        }
      ]
    }
  ]
}
```

## 📄 **TSV Receipt Example**

The TSV format provides a row-by-row breakdown of each process for every analysis, which makes it easy to filter and analyze in spreadsheet applications:

```tsv
submitter_analysis_id	overall_status	process	status	exit_code	timestamp	error_message	file_manager_analysis_id	analysis_type	study_id	analysis_state	published_at
analysis_004	FAILED	CHECK_SUBMISSION_DEPENDENCIES:ANALYSIS_SPLIT	PASS	0	2025-08-28 20:51:01.287276		Not applicable	sequenceAlignment	TEST-CA	Not applicable	Not applicable
analysis_004	FAILED	PCGL:MOLECULAR_DATA_SUBMISSION_WORKFLOW:CHECK_SUBMISSION_DEPENDENCIES:VALIDATE_CLINICAL	SUCCESS	0	2025-08-28T20:51:09+00:00		Not applicable	sequenceAlignment	TEST-CA	Not applicable	Not applicable
analysis_004	FAILED	PCGL:MOLECULAR_DATA_SUBMISSION_WORKFLOW:CHECK_SUBMISSION_DEPENDENCIES:CLINICAL_SUBMISSION	SUCCESS	0	2025-08-28T20:51:14+00:00		Not applicable	sequenceAlignment	TEST-CA	Not applicable	Not applicable
analysis_004	FAILED	PCGL:MOLECULAR_DATA_SUBMISSION_WORKFLOW:METADATA_PAYLOAD_GENERATION:PAYLOAD_GENERATE	FAILED	1	2025-08-28T20:51:18+00:00	Error: workflow-meta is required for analysis type 'sequenceAlignment'	Not applicable	sequenceAlignment	TEST-CA	Not applicable	Not applicable
analysis_005	SUCCESS	CHECK_SUBMISSION_DEPENDENCIES:ANALYSIS_SPLIT	PASS	0	2025-08-28 20:52:01.123456		0fe74e99-f30e-40b3-a74e-99f30eb0b302	variantCall	TEST-CA	PUBLISHED	2025-08-28T20:55:00Z
analysis_005	SUCCESS	PCGL:MOLECULAR_DATA_SUBMISSION_WORKFLOW:METADATA_PAYLOAD_GENERATION:PAYLOAD_GENERATE	SUCCESS	0	2025-08-28T20:52:15+00:00		0fe74e99-f30e-40b3-a74e-99f30eb0b302	variantCall	TEST-CA	PUBLISHED	2025-08-28T20:55:00Z
analysis_005	SUCCESS	PCGL:MOLECULAR_DATA_SUBMISSION_WORKFLOW:SONG_SUBMISSION:SONG_SUBMIT	SUCCESS	0	2025-08-28T20:54:20+00:00		0fe74e99-f30e-40b3-a74e-99f30eb0b302	variantCall	TEST-CA	PUBLISHED	2025-08-28T20:55:00Z
analysis_005	SUCCESS	PCGL:MOLECULAR_DATA_SUBMISSION_WORKFLOW:SONG_SUBMISSION:SONG_PUBLISH	SUCCESS	0	2025-08-28T20:54:58+00:00		0fe74e99-f30e-40b3-a74e-99f30eb0b302	variantCall	TEST-CA	PUBLISHED	2025-08-28T20:55:00Z
```

### **TSV Format Benefits:**
- **One row per process**: Each workflow process gets its own row for detailed tracking
- **Easy filtering**: Filter by analysis ID, process name, or status in spreadsheet applications  
- **Error tracking**: Quickly identify which specific processes failed for each analysis

## 🔍 **Understanding Process Status Values**

| Status | Meaning | Action Required |
|--------|---------|-----------------|
| `SUCCESS` | Process completed successfully | None - continue to next step |
| `FAILED` | Process encountered an error | Review error_message and correct input data |

## 🛠️ **Troubleshooting with Receipts**

### **Common Error Patterns**

When reviewing failed processes in your batch receipt, you'll encounter various error patterns. For detailed descriptions of common error types and their solutions, see the [Troubleshooting Guide (TBD)](troubleshoot.md).

### **Using Work Directories for Debugging**

Each process includes a `work_directory` field that points to the execution directory of the process. This directory contains:
- Input files used by the process
- Output files generated from the process
- Log files with detailed execution information

### **How to Resolve Issues**

1. **Identify Failed Processes**: Look for `"status": "FAILED"` in the JSON receipt
2. **Read Error Messages**: Check the `error_message` field in process details
3. **Access Work Directory**: Navigate to the `work_directory` path for detailed logs
4. **Correct Input Data**: Fix the identified issues in your metadata or data files
5. **Re-run Workflow**: Submit the corrected data through the workflow again

## 💡 **Best Practices**

- Always review the batch receipt after workflow completion
- Save receipts for audit purposes as they contain complete processing history
- Choose either JSON or TSV format receipt according to your own preference
- Check work directories when error messages need more context
