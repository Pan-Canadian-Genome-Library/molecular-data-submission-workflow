# PCGL Molecular Data Submission Workflow: Usage Guide

This guide provides practical instructions and command-line examples for running the PCGL molecular data submission workflow.

## 🚀 Quick Start

### **Prerequisites Checklist**
Before running the workflow, ensure you have:
- [ ] Study registered in PCGL (Contact: helpdesk@genomelibrary.ca)
- [ ] All participants registered (Coordinate with your Study data coordinator)
- [ ] API token obtained (Request from PCGL administrator)
- [ ] Prepared metadata files (see **[Input Documentation](Input.md)** for requirements and formats)
- [ ] Data files organized in accessible directory
- [ ] Nextflow installed (version 23.04.1 or later)
- [ ] Container engine (Docker recommended)

### **Minimal Required Command**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "STUDY_001" \
    --token "your_token_here" \
    --path_to_files_directory "/path/to/data/files" \
    --file_metadata "metadata/files.tsv" \
    --analysis_metadata "metadata/analyses.tsv" \
    --outdir results \
    -profile docker,sd4h_prod
```

## 🏭 Production Usage Examples

### **1. Complete Workflow (All Metadata)**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "EXAMPLE-CA9" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/data/genomics/study_files" \
    --file_metadata "metadata/file_metadata.tsv" \
    --analysis_metadata "metadata/analysis_metadata.tsv" \
    --workflow_metadata "metadata/workflow_metadata.tsv" \
    --read_group_metadata "metadata/read_group_metadata.tsv" \
    --experiment_metadata "metadata/experiment_metadata.tsv" \
    --specimen_metadata "metadata/specimen_metadata.tsv" \
    --sample_metadata "metadata/sample_metadata.tsv" \
    --outdir submission_results \
    -profile docker,sd4h_prod
```

### **2. Minimal Workflow (Required Files Only)**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "EXAMPLE-CA9" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/data/genomics/study_files" \
    --file_metadata "metadata/file_metadata.tsv" \
    --analysis_metadata "metadata/analysis_metadata.tsv" \
    --outdir submission_results \
    -profile docker,sd4h_prod
```

### **3. Partial Metadata Submission**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "EXAMPLE-CA9" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/data/genomics/study_files" \
    --file_metadata "metadata/file_metadata.tsv" \
    --analysis_metadata "metadata/analysis_metadata.tsv" \
    --experiment_metadata "metadata/experiment_metadata.tsv" \
    --sample_metadata "metadata/sample_metadata.tsv" \
    --outdir submission_results \
    -profile docker,sd4h_prod
    # Note: workflow_metadata, read_group_metadata, specimen_metadata omitted
```

### **4. Debug Mode with Detailed Logging**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "EXAMPLE-CA9" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/data/genomics/study_files" \
    --file_metadata "metadata/file_metadata.tsv" \
    --analysis_metadata "metadata/analysis_metadata.tsv" \
    --debug_channels true \
    --outdir debug_results \
    -profile docker \
    -with-trace \
    -with-timeline \
    -with-report
```

## ⚙️ Environment Configuration

### **PCGL Environment Profiles**

The workflow supports different PCGL environments through configuration profiles:

#### **Development Environment**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile docker,cumulus_dev
```

#### **QA Environment**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile docker,sd4h_qa
```

#### **Production Environment**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile docker,sd4h_prod
```

### **Container Engine Options**

#### **Docker (Recommended)**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile docker
```

#### **Singularity**
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile singularity
```


## 📚 Related Documentation

- **[Input Documentation](Input.md)** - Parameter requirements, metadata schemas, and file formats
- **[Testing Documentation](testing.md)** - Testing procedures with sample data
- **[Output Documentation](output.md)** - Understanding workflow results and receipts
- **[Troubleshooting Documentation](troubleshooting.md)** - Common issues and solutions
- **[Receipt Documentation](receipt.md)** - Interpreting batch receipts and status information