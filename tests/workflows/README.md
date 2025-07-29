# Molecular Data Submission Workflow Tests

This directory contains comprehensive test cases for the `MOLECULAR_DATA_SUBMISSION_WORKFLOW`.

## Test Structure

### Test Data (`test_data/`)
- **Metadata files**: CSV files containing mock metadata for different entities
  - `file_metadata.csv`: File information and checksums
  - `analysis_metadata.csv`: Analysis details and types
  - `workflow_metadata.csv`: Workflow execution metadata
  - `read_group_metadata.csv`: Sequencing read group information
  - `experiment_metadata.csv`: Experimental design details
  - `specimen_metadata.csv`: Biological specimen information
  - `sample_metadata.csv`: Sample and donor mapping

- **Data files** (`files/`): Mock data files for testing
  - `test_file_001_1.bam`, `test_file_001_2.bam`: BAM files for analysis_001
  - `test_file_002_1.bam`: BAM file for analysis_002  
  - `test_file_003_1.vcf`: VCF file for analysis_003

### Test Cases

#### 1. Complete Workflow Test (`Should run complete workflow with upload enabled`)
- Tests the full workflow with all subworkflows enabled
- Upload is enabled (skip_upload = false)
- Duplicate check is disabled for faster testing
- Validates successful completion and output channels

#### 2. Upload Skipped Test (`Should run workflow with upload skipped`)
- Tests workflow behavior when data upload is skipped
- Key difference: skip_upload = true
- Ensures workflow handles empty analysis channel gracefully
- Tests receipt generation with null analysis JSON

#### 3. Duplicate Check Test (`Should handle duplicate check enabled`)
- Tests workflow with duplicate checking enabled
- Validates that duplicate detection doesn't break the pipeline
- Ensures all processing steps complete successfully

#### 4. Analysis Types Test (`Should process different analysis types`)
- Tests handling of different genomics analysis types:
  - `sequenceExperiment`: Raw sequencing data
  - `sequenceAlignment`: Aligned sequence data  
  - `variantCall`: Variant calling results
- Validates type-specific processing and file handling

#### 5. Minimal Configuration Test (`Should run minimal workflow`)
- Tests workflow with minimal configuration
- Both duplicate check and upload are skipped
- Validates core workflow functionality without optional components
- Useful for CI/CD pipelines and quick validation

## Test Data Design

The test data includes:
- **3 analyses** with different types and statuses
- **Mixed success/failure scenarios** to test error handling
- **Realistic file types** (BAM, VCF) for different analysis types
- **Complete metadata chain** from files to samples to donors

## Running Tests

```bash
# Run all tests
nf-test test tests/workflows/molecular-data-submission-workflow.nf.test

# Run specific test
nf-test test tests/workflows/molecular-data-submission-workflow.nf.test -t "Should run complete workflow with upload enabled"

# Run with verbose output
nf-test test tests/workflows/molecular-data-submission-workflow.nf.test -v
```

## Expected Behavior

All tests should:
1. ✅ Complete successfully (`workflow.success`)
2. ✅ Generate version information (`workflow.out.versions`)
3. ✅ Create status files for tracking (`workflow.out.all_status`)
4. ✅ Execute multiple tasks (`workflow.trace.tasks().size() > 0`)

## Test Coverage

These tests cover:
- 🧪 **Core workflow logic**: All major subworkflows and processing steps
- 🧪 **Configuration variants**: Different parameter combinations
- 🧪 **Error scenarios**: Failed analyses and status propagation
- 🧪 **Data types**: Multiple genomics file formats and analysis types
- 🧪 **Channel operations**: Grouping, joining, and filtering logic
- 🧪 **Receipt generation**: Batch receipt creation and aggregation
