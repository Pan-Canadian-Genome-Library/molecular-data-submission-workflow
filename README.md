# Pan-Canadian-Genome-Library/molecular-data-submission-workflow

## Introduction

**Pan-Canadian-Genome-Library/molecular-data-submission-workflow** is a Nextflow pipeline that automates the validation, packaging, and submission of molecular genomics data and associated metadata to the Pan-Canadian Genome Library (PCGL) data repository. The pipeline ensures data integrity, validates metadata compliance, handles file uploads, and generates comprehensive submission receipts for tracking and audit purposes.

## Pipeline Overview

The workflow consists of five main stages:

1. **Dependency Checking** - Validates input files, metadata completeness, and submission prerequisites
2. **Metadata Payload Generation** - Creates standardized JSON payloads from input metadata
3. **Data Validation** - Validates data files and metadata against PCGL schemas and requirements  
4. **Data Upload** - Handles secure file transfer and metadata submission to PCGL repositories
5. **Receipt Generation** - Creates detailed batch receipts and summary reports for submission tracking

## Prerequisites

### System Requirements

- **Nextflow**: Version 22.04.0 or newer with DSL2 support
- **Container Engine**: One of the following:
  - Docker (recommended)
  - Singularity/Apptainer 
  - Conda (alternative, but containers preferred)
- **Java**: Java 8 or later (required by Nextflow)
- **Memory**: Minimum 8GB RAM recommended
- **Storage**: Sufficient disk space for input data, intermediate files, and outputs

### Access Requirements

- **PCGL Access Token**: Valid authentication token with submission permissions for your study
- **Study Registration**: Your study must be registered in the PCGL system
- **Participant Registration**: The participants in your study must be registered in the PCGL system
- **Network Access**: Connectivity to PCGL submission endpoints:
  - File Manager service
  - File Transfer service  
  - Clinical submission service 

### Input Data Requirements

- **Molecular Data Files**: 
  - Supported formats: CRAM, BAM, VCF, BCF
  - Files must include appropriate index files (e.g., .crai, .bai, .tbi)
  - Files must be accessible from the specified `path_to_files_directory`

- **Metadata Files**: 
  - **Required**: `file_metadata.tsv`, `analysis_metadata.tsv`
  - **Optional**: `workflow_metadata.tsv`, `read_group_metadata.tsv`, `experiment_metadata.tsv`, `specimen_metadata.tsv`, `sample_metadata.tsv`
  - All metadata files must be in tab-separated (TSV) format
  - Files must comply with the your study cumstom data model, which is the combination of PCGL Base and Extentions Data Model. For the latest version of the PCGL Base Data Model, please see the [latest release folder](https://drive.google.com/drive/u/1/folders/1vfNA7ajwh3WKkbVmswb6j9TuWKxaN9bB).

### Environment Setup

1. **Install Nextflow**:
   ```bash
   curl -s https://get.nextflow.io | bash
   mv nextflow ~/bin/ # or add to your PATH
   ```

2. **Install Container Engine**:
   - **Docker**: Follow [Docker installation guide](https://docs.docker.com/get-docker/)
   - **Singularity**: Follow [Singularity installation guide](https://sylabs.io/guides/3.0/user-guide/installation.html)

3. **Verify Installation**:
   ```bash
   nextflow info
   docker --version  # or singularity --version
   ```

4. **Obtain PCGL Access Token**:
   - Contact your PCGL administrator or data coordinator
   - Ensure the token has appropriate permissions for your study

## Usage

> [!NOTE]
> If you are new to Nextflow, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to test your setup before running the workflow on actual data.

First, prepare your metadata files and data directory structure, e.g:

```
input/
├── metadata/
│   ├── file_metadata.tsv          # Required: File information  
│   ├── analysis_metadata.tsv      # Required: Analysis details
│   ├── workflow_metadata.tsv      # Optional: Workflow information
│   ├── read_group_metadata.tsv    # Optional: Read group details
│   ├── experiment_metadata.tsv    # Optional: Experiment information
│   ├── specimen_metadata.tsv      # Optional: Specimen details
│   └── sample_metadata.tsv        # Optional: Sample information
└── data/
    └── [your data files - CRAM, BAM, etc.]
```

### Basic Usage with Required Metadata

```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "YOUR_STUDY_ID" \
    --token "YOUR_ACCESS_TOKEN" \
    --path_to_files_directory "/path/to/data" \
    --file_metadata "/path/to/file_metadata.tsv" \
    --analysis_metadata "/path/to/analysis_metadata.tsv" \
    --outdir results \
    -profile docker
```

### Advanced Usage with Additional Metadata

```bash
nextflow run Pan-Canadian-Genome-Library/molecular-data-submission-workflow \
    --study_id "YOUR_STUDY_ID" \
    --token "YOUR_ACCESS_TOKEN" \
    --path_to_files_directory "/path/to/data" \
    --file_metadata "metadata/file_metadata.tsv" \
    --analysis_metadata "metadata/analysis_metadata.tsv" \
    --workflow_metadata "metadata/workflow_metadata.tsv" \
    --read_group_metadata "metadata/read_group_metadata.tsv" \
    --experiment_metadata "metadata/experiment_metadata.tsv" \
    --specimen_metadata "metadata/specimen_metadata.tsv" \
    --sample_metadata "metadata/sample_metadata.tsv" \
    --outdir results \
    -profile docker
```

For more detailed usage instructions, see [docs/usage.md](docs/usage.md).

## Testing

Before running the pipeline on your data, we recommend testing it with the provided test datasets:

- **[Testing Guide](docs/Testing.md)** - Comprehensive instructions on how to test the workflow using `nextflow run` with the included test datasets
- Test datasets are provided in the `tests/test_data/` directory
- Both minimal and comprehensive testing scenarios are covered


## Credits

Pan-Canadian-Genome-Library/molecular-data-submission-workflow was originally written by edsu7,lxiang.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use Pan-Canadian-Genome-Library/molecular-data-submission-workflow for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
