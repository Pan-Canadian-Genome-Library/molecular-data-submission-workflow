# PCGL Molecular Data Submission Workflow: Introduction

## Overview

The Pan-Canadian Genome Library (PCGL) Molecular Data Submission Workflow is a Nextflow-based pipeline designed to automate the validation, packaging, and submission of molecular genomics data and associated metadata to the PCGL data repository. The workflow ensures data integrity, metadata compliance, secure file transfer, and comprehensive submission tracking.

This document provides a detailed overview of the workflow components, prerequisites, and system architecture to help users understand, configure, and operate the pipeline effectively.

---

## Workflow Components

The workflow is organized into five main subworkflows, each responsible for a critical stage of the submission process:

### 1. Dependency Checking
- Validates study and participant registration in PCGL systems
- Assesses biospecimen entity status (samples, specimens, experiments, read groups)
- Verifies existence and accessibility of required input files
- Cross-validates metadata consistency with existing PCGL records
- Registers missing biospecimen entities as needed

### 2. Metadata Payload Generation
- Converts input metadata from TSV to standardized JSON payloads
- Validates payloads against PCGL data model JSON schemas

### 3. Data Validation
- Performs format and completeness checks on metadata files
- Assesses sequencing data integrity and format compliance
- Cross-validates relationships between sequencing data and metadata

### 4. Data Uploading
- Submits validated JSON payloads to PCGL file manager
- Generates file manifest with checksums and transfer metadata
- Uploads sequencing files to PCGL object storage
- Publishes analysis metadata and retrieves unique analysis IDs

### 5. Receipt Generation
- Provides real-time submission status feedback
- Generates comprehensive batch receipts in JSON and TSV formats
- Summarizes results, errors, and audit information

---

## Prerequisites

### System Requirements
- **Nextflow**: v22.04.0+ (DSL2)
- **Container Engine**: Docker (recommended), Singularity/Apptainer, or Conda
- **Java**: 17 or newer
- **Bash**: 3.2 or newer
- **Memory**: 8GB RAM minimum
- **Storage**: Sufficient disk space for input, intermediate, and output files

### Access Requirements
- **PCGL API Token**: Valid token with submission permissions
- **Network Access**: Connectivity to PCGL File Manager, File Transfer, and Clinical Submission services

### Submission Dependencies
- **Study Registration**: Study must be registered in PCGL
- **Participant Registration**: All participants in your submission batch must be registered
- **Biospecimen Entities**: Metadata for samples, specimens, experiments, and read groups must be provided if not already submitted

### Input Data Requirements
- **Molecular Data Files**: CRAM, BAM, VCF, BCF (with index files)
- **Metadata Files**: Required: `file_metadata.tsv`, `analysis_metadata.tsv`; Optional: biospecimen and workflow metadata files
- **Format Compliance**: All metadata must follow PCGL Base and Extension Data Model specifications

---

## System Architecture

The workflow leverages Nextflow DSL2 for modular, reproducible, and scalable execution. Key architectural features include:

- **Modular Subworkflows**: Each stage is implemented as a separate module for maintainability and extensibility
- **Containerization**: All processes run in containers for reproducibility and portability
- **Centralized Configuration**: Pipeline parameters and environment settings are managed via configuration files
- **Integrated Validation**: Metadata and data validation are performed at multiple stages to ensure compliance
- **Secure Data Transfer**: Genomic files and metadata are uploaded using secure protocols to PCGL endpoints
- **Comprehensive Logging**: All steps generate logs and status files for audit and troubleshooting
- **Receipt-Based Tracking**: Batch receipts provide hierarchical tracking of submission status at batch, analysis, and process levels

---

## References & Further Reading
- [PCGL Base Data Model](https://drive.google.com/drive/u/1/folders/1vfNA7ajwh3WKkbVmswb6j9TuWKxaN9bB)
- [nf-core Framework](https://nf-co.re/)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)
- [PCGL Molecular Data Submission Workflow GitHub](https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow)

---

For detailed usage instructions, troubleshooting, and output interpretation, see the corresponding documentation files:
- [Input Documentation](input.md)
- [Testing Guide](testing.md)
- [Receipt Guide](receipt.md)
- [Output Documentation](output.md)
- [Troubleshooting Guide](troubleshooting.md)

---

For questions or support, please contact your PCGL administrator or data coordinator.