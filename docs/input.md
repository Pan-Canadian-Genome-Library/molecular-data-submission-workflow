# Input

## Required Parameters

| Parameter                  | Description                                           | Type     |
|----------------------------|-------------------------------------------------------|----------|
| `study_id`                 | PCGL study identifier                                 | string   |
| `token`                    | Authentication token with PCGL submission scope      | string   |
| `path_to_files_directory`  | Path to directory containing files to be uploaded    | path     |
| `file_metadata`            | Tab-separated file containing file information        | file     |
| `analysis_metadata`        | Tab-separated file containing analysis information    | file     |

## Optional Parameters

| Parameter                  | Description                                           | Default  |
|----------------------------|-------------------------------------------------------|----------|
| `workflow_metadata`        | Tab-separated file containing workflow information    | null     |
| `read_group_metadata`      | Tab-separated file containing read group information  | null     |
| `experiment_metadata`      | Tab-separated file containing experiment information  | null     |
| `specimen_metadata`        | Tab-separated file containing specimen information    | null     |
| `sample_metadata`          | Tab-separated file containing sample information      | null     |

## Advanced Configuration Parameters

| Parameter                           | Description                                      | Default           |
|-------------------------------------|--------------------------------------------------|-------------------|
| `outdir`                            | Output directory for results                     | "out"             |
| `skip_upload`                       | Skip actual molecular data upload.               | false             |
| `allow_duplicates`                  | Allow duplicate file submissions                 | false             |
| `exit_on_error`                     | Exit pipeline on first error                    | false             |
| `debug_channels`                    | Enable debug output for workflow channels       | false             |
| `file_manager_url`                  | File manager server URL                   | null              |
| `file_transfer_url`                 | File transfer server URL                  | null              |
| `clinical_url`                      | Clinical submission URL                   | null              |
| `file_manager_container`            | File manager client container image                     | ghcr.io/overture-stack/song-client |
| `file_manager_container_tag`        | File manager client container tag                       | 822055be          |
| `file_transfer_container`           | File transfer container image                    | ghcr.io/pan-canadian-genome-library/file-transfer |
| `file_transfer_container_tag`       | File transfer container tag                     | edge              |
| `file_transfer_transport_parallel`  | File transfer parallel configuration                  | null              |
| `file_transfer_transport_mem`       | File transfer memory configuration               | null              |

## Input File Validation

The workflow automatically validates:
- File existence and accessibility
- Metadata file format and required columns
- Data file integrity and format compliance
- Metadata consistency and relationships

## File Format Specifications

### Molecular Data Files
- **Supported formats**: CRAM, BAM, VCF, BCF
- **Index files**: Required for each data file (e.g., .crai, .bai, .tbi, .csi)
- **Naming**: File names must match entries in `file_metadata.tsv`
- **Location**: All files must be accessible from the specified `path_to_files_directory`
- **Integrity**: Files should be properly formatted and not corrupted

### Metadata File Schema Examples

#### `file_metadata.tsv` (Required)
| Column                | Description                        | Example                | Required |
|-----------------------|------------------------------------|------------------------|----------|
| submitter_analysis_id | Unique analysis identifier         | analysis_001           | Yes      |
| fileName              | Name of the data file              | sample001.bam          | Yes      |
| fileSize              | Size in bytes                      | 123456789              | Yes      |
| fileMd5sum            | MD5 checksum of the file           | 1a2b3c4d5e6f...        | Yes      |
| fileType              | Type of file (CRAM/BAM/VCF/BCF)    | BAM                    | Yes      |
| fileAccess            | File access level                  | controlled             | Yes      |
| dataType              | Type of data content               | Aligned Reads          | Yes      |

#### `analysis_metadata.tsv` (Required)
| Column                     | Description                        | Example                | Required |
|----------------------------|------------------------------------|------------------------|----------|
| studyId                    | Study identifier                   | TEST-CA                | Yes      |
| submitter_analysis_id      | Unique analysis identifier         | analysis_001           | Yes      |
| analysisType               | Type of analysis performed         | sequenceAlignment      | Yes      |
| submitter_participant_id   | Participant identifier             | PART_001               | No      |
| submitter_specimen_id      | Specimen identifier                | SPEC_001               | No       |
| submitter_sample_id        | Sample identifier                  | SAMP_001               | No       |
| submitter_experiment_id    | Experiment identifier              | EXP_001                | No       |
| data_category              | Category of data                   | Sequencing Reads       | Yes      |
| variant_class              | Class of variants                  | SNV                    | No       |
| variant_calling_strategy   | Strategy used for variant calling  | WGS                    | No       |
| genome_build               | Reference genome build             | GRCh38                 | No       |
| genome_annotation          | Genome annotation version          | GENCODE v29            | No       |

#### `workflow_metadata.tsv` (Optional)
| Column                 | Description                        | Example                | Required |
|------------------------|------------------------------------|------------------------|----------|
| submitter_workflow_id  | Unique workflow identifier         | workflow_001           | Yes      |
| submitter_analysis_id  | Linked analysis identifier         | analysis_001           | Yes      |
| workflow_name          | Name of workflow                   | bwa-mem2-alignment     | Yes      |
| workflow_version       | Version of workflow                | 2.2.1                  | No      |
| workflow_url           | URL to workflow repository         | github.com/...         | No       |

#### Biospecimen Metadata Files (Optional)
For biospecimen metadata files (`specimen_metadata.tsv`, `sample_metadata.tsv`, `experiment_metadata.tsv`, `read_group_metadata.tsv`):
- Each file must include unique entity IDs as primary keys
- All required fields must be present as per PCGL Base Data Model
- Foreign key relationships must be maintained between entities
- See [PCGL Base Data Model](https://drive.google.com/drive/u/1/folders/1vfNA7ajwh3WKkbVmswb6j9TuWKxaN9bB) for complete schemas

## Data Preparation Guidelines

- Validate all metadata files for required columns and correct formats
- Ensure data files are properly indexed and not corrupted
- Use consistent naming conventions across metadata and data files
- Verify file checksums match actual file content
- Ensure data files are accessible from the specified `path_to_files_directory`
- Test with small batch before large submissions

### Recommended Directory Structure

This suggested organization separates metadata files from data files, making it easier to manage and reference your submission components:

```
study/
├── metadata/
│   ├── file_metadata.tsv          # Required: File information
│   ├── analysis_metadata.tsv      # Required: Analysis information
│   ├── workflow_metadata.tsv      # Optional: Workflow information
│   ├── read_group_metadata.tsv    # Optional: Read group information
│   ├── experiment_metadata.tsv    # Optional: Experiment information
│   ├── specimen_metadata.tsv      # Optional: Specimen information
│   └── sample_metadata.tsv        # Optional: Sample information
└── data/                          # path_to_files_directory
    ├── sample1.cram
    ├── sample1.cram.crai
    ├── sample2.cram
    └── sample2.cram.crai
```

**Benefits of this structure:**
- **Clear separation**: Metadata and data files are organized in separate directories
- **Easy maintenance**: All TSV metadata files are centralized in one location
- **Simple referencing**: The `data/` directory serves as your `path_to_files_directory` parameter
- **Scalability**: Easy to add more files without cluttering the project root
- **Index file co-location**: Data files and their corresponding index files are kept together

---

For further details, see:
- [Introduction Guide](Introduction.md) - Workflow overview and architecture
- [Testing Guide](Testing.md) - Testing instructions with sample data
- [Receipt Guide](receipt.md) - Understanding submission results
- [Output Documentation](output.md) - Complete output file descriptions
- [Troubleshooting Guide](troubleshooting.md) - Common issues and solutions