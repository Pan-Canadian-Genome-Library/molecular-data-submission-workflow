# PCGL Molecular Data Submission Workflow: Usage Guide

This guide provides practical instructions and command-line examples for running the PCGL molecular data submission workflow.

## 🚀 Quick Start

### **Prerequisites Checklist**
Before running the workflow, ensure you have:
- [ ] Study registered in PCGL (Contact: helpdesk@genomelibrary.ca)
- [ ] All participants registered (Coordinate with your Study data coordinator)
- [ ] API token obtained (Contact: helpdesk@genomelibrary.ca)
- [ ] Prepared metadata files (see **[Input Documentation](Input.md)** for requirements and formats)
- [ ] Data files organized in accessible directory
- [ ] Nextflow installed (version 22.04.2 or later)
- [ ] Container engine installed (Docker or Singularity)

### **Minimal Required Command**
Submission with only the essential metadata files required to start a submission:
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "your_pcgl_study_id" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/path/to/data/files" \
    --file_metadata "/path/to/file_metadata" \
    --analysis_metadata "/path/to/analysis_metadata" \
    --outdir "/path/to/submission_results" \
    -profile [docker|singularity],sd4h_prod
```

## 🏭 Production Usage Examples

### **1. All Metadata Submission**
Submission including all optional metadata files provided:
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "your_pcgl_study_id" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/path/to/data/files" \
    --file_metadata "/path/to/file_metadata" \
    --analysis_metadata "/path/to/analysis_metadata" \
    --workflow_metadata "/path/to/workflow_metadata" \
    --read_group_metadata "/path/to/read_group_metadata" \
    --experiment_metadata "/path/to/experiment_metadata" \
    --specimen_metadata "/path/to/specimen_metadata" \
    --sample_metadata "/path/to/sample_metadata" \
    --outdir "/path/to/submission_results" \
    -profile [docker|singularity],sd4h_prod
```

### **2. Minimal Metadata Submission**
Submission with only the two required metadata files (files and analyses):
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "your_pcgl_study_id" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/path/to/data/files" \
    --file_metadata "/path/to/file_metadata" \
    --analysis_metadata "/path/to/analysis_metadata" \
    --outdir "/path/to/submission_results" \
    -profile [docker|singularity],sd4h_prod
```

### **3. Partial Metadata Submission**
Submission with selective metadata files - include only what's available or relevant:
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "your_pcgl_study_id" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/path/to/data/files" \
    --file_metadata "/path/to/file_metadata" \
    --analysis_metadata "/path/to/analysis_metadata" \
    --experiment_metadata "/path/to/experiment_metadata" \
    --sample_metadata "/path/to/sample_metadata" \
    --outdir "/path/to/submission_results" \
    -profile [docker|singularity],sd4h_prod
    # Note: workflow_metadata, read_group_metadata, specimen_metadata omitted
```

### **4. Enable Detailed Channel Output Logging**
Enable verbose logging to troubleshoot issues or monitor data flow between channels:
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "your_pcgl_study_id" \
    --token "your_pcgl_token_here" \
    --path_to_files_directory "/path/to/data/files" \
    --file_metadata "/path/to/file_metadata" \
    --analysis_metadata "/path/to/analysis_metadata" \
    --debug_channels true \
    --outdir "/path/to/submission_results" \
    -profile [docker|singularity],sd4h_prod 
```

## ⚙️ Environment Configuration

### **PCGL Environment Profiles**

The workflow supports different PCGL environments through configuration profiles:

#### **PCGL SD4H Development Environment**
Connect to the development environment for testing and validation:
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile sd4h_dev
```

#### **PCGL SD4H QA Environment**
Connect to the quality assurance environment for pre-production testing:
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile sd4h_qa
```

#### **PCGL SD4H Production Environment**
Connect to the production environment for live data submissions:
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile sd4h_prod
```

### **Container Engine Options**

#### **Docker**
Use Docker containers (recommended for most systems with Docker installed):
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile docker
```

#### **Singularity**
Use Singularity containers (recommended for HPC environments without Docker):
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile singularity
```

### **Compute Cluster Management Options**

#### **Slurm**
Use Slurm  (recommended for HPC environments where available and required for resource management):
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    -profile slurm
```

Note the config provided in :
```
conf/slurm.config
```
Is rudimentary and serves as a start/template, likely to not work out of box. Please add in configurations and arguements required (such as specific nodes).


## 📚 Related Documentation

- **[Input Documentation](Input.md)** - Parameter requirements, metadata schemas, and file formats
- **[Testing Documentation](testing.md)** - Testing procedures with sample data
- **[Output Documentation](output.md)** - Understanding workflow results and receipts
- **[Troubleshooting Documentation](troubleshooting.md)** - Common issues and solutions
- **[Receipt Documentation](receipt.md)** - Interpreting batch receipts and status information