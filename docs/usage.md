# PCGL Molecular Data Submission Workflow: Usage Guide

This guide provides practical instructions and command-line examples for running the PCGL molecular data submission workflow.

## 🚀 Quick Start

### **Prerequisites Checklist**
Before running the workflow, ensure the following are in place:
- [ ] **Study registered** — Your study is registered in PCGL (contact [helpdesk@genomelibrary.ca](mailto:helpdesk@genomelibrary.ca))
- [ ] **Participants registered** — All participants in the submission batch are registered (coordinate with your Study Data Coordinator)
- [ ] **API token obtained** — A valid API token has been issued for your study (contact [helpdesk@genomelibrary.ca](mailto:helpdesk@genomelibrary.ca))
- [ ] **Metadata files prepared** — All required metadata files conform to the data model (see **[Input Documentation](input.md)** for schemas and formatting requirements)
- [ ] **Data files accessible** — Molecular data files are organized in a readable directory or referenced by absolute paths
- [ ] **Nextflow installed** — Version 22.04.2 or later ([installation guide](https://www.nextflow.io/docs/latest/install.html))
- [ ] **Container engine installed** — Docker or Singularity is available on the target system

### **Minimal Required Command**
Submission with only the essential metadata files required to start a submission (`path_to_files_directory` is optional when `fileName` in `file_metadata.tsv` contains absolute paths):
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

> **File pathing**: `path_to_files_directory` can be omitted when the `fileName` column in `file_metadata.tsv` contains absolute paths or subdirectory-relative paths that already resolve to accessible files. See [Input Documentation](input.md#molecular-data-files) for the full set of supported pathing options.

> **Checksums**: `fileSize` and `fileMd5sum` columns in `file_metadata.tsv` are **optional**. The workflow calculates them automatically from the actual files and embeds the values in the submission payload.

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

### **5. Submission Using Absolute File Paths (no `path_to_files_directory`)**
When data files are stored at absolute paths (e.g., on a shared storage mount), you can omit `--path_to_files_directory` and reference them directly in `file_metadata.tsv`:
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "your_pcgl_study_id" \
    --token "your_pcgl_token_here" \
    --file_metadata "/path/to/file_metadata" \
    --analysis_metadata "/path/to/analysis_metadata" \
    --outdir "/path/to/submission_results" \
    -profile [docker|singularity],sd4h_prod
    # Note: fileName in file_metadata.tsv must contain absolute paths, e.g. /mnt/storage/sample001.cram
```

#### **Adjusting max forks**
Resource management is important when submitting a batch of files. By default we conservatitely limit the number of `PAYLOAD_GENERATE` and `SCORE_UPLOAD` tasks to the default value of `1`. This means only 1 of each task can be running at a time.\n
This is to reduce I/O overhead and ensure efficient transfer.  If the data is local and not constrained, the number can increased.
```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    [... parameters ...] \
    --fork_limit 1
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
Is rudimentary and serves as a start/template, likely to not work out of box. Please add in configurations and arguments required (such as specific nodes).


## 📚 Related Documentation

- **[Input Documentation](input.md)** - Parameter requirements, metadata schemas, and file formats
- **[Testing Documentation](testing.md)** - Testing procedures with sample data
- **[Output Documentation](output.md)** - Understanding workflow results and receipts
- **[Troubleshooting Documentation (TBD)](troubleshooting.md)** - Common issues and solutions
- **[Receipt Documentation](receipt.md)** - Interpreting batch receipts and status information