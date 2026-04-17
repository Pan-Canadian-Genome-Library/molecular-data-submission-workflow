# PCGL Molecular Data Submission Workflow: Input

## Quick Start Checklist

Before running the workflow, ensure you have:
- [ ] `study_id` and `token` from PCGL administrator  
- [ ] Completed `file_metadata.tsv` and `analysis_metadata.tsv`
- [ ] Data files with corresponding index files in accessible directory
- [ ] Verified file names in metadata match actual file names 

**New to PCGL?** See [Getting Started](introduction.md#getting-started) for registration steps.  
**Ready to run?** See [Usage Guide](usage.md) for execution commands.

## Required Parameters

| Parameter                  | Description                                           | Type     |
|----------------------------|-------------------------------------------------------|----------|
| `study_id`                 | PCGL study identifier                                 | string   |
| `token`                    | Authentication token with PCGL submission scope      | string   |
| `file_metadata`            | Tab-separated file containing file information        | file     |
| `analysis_metadata`        | Tab-separated file containing analysis information    | file     |



## Optional Parameters

| Parameter                  | Description                                           | Type  |
|----------------------------|-------------------------------------------------------|----------|
| `workflow_metadata`        | Tab-separated file containing workflow information    | file     |
| `read_group_metadata`      | Tab-separated file containing read group information  | file     |
| `experiment_metadata`      | Tab-separated file containing experiment information  | file     |
| `specimen_metadata`        | Tab-separated file containing specimen information    | file     |
| `sample_metadata`          | Tab-separated file containing sample information      | file     |
| `path_to_files_directory`  | Path to directory containing files to be uploaded    | string     |

> **Note on file pathing**: `path_to_files_directory` is optional. When omitted, the `fileName` column in `file_metadata.tsv` must contain the absolute path to each file.

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

## File Requirements and Specifications

### Molecular Data Files
- **Supported formats**: CRAM, BAM, VCF, BCF
- **Index files**: Required for each data file (e.g., .crai, .bai, .tbi, .csi)
- **Naming / path requirements**:
  - The `fileName` column in `file_metadata.tsv` accepts:
    - A plain file name (e.g., `sample001.bam`) — file is looked up relative to `path_to_files_directory`
    - A subdirectory path (e.g., `wgs/sample001.bam`) — resolved under `path_to_files_directory`
    - An absolute path (e.g., `/data/submissions/sample001.bam`) — used directly, `path_to_files_directory` is ignored for that file
  - Regardless of the path format used, **only the base filename is retained** by the workflow (directory components are stripped). For example, both `lane1/sample.bam` and `/data/lane2/sample.bam` resolve to `sample.bam`. Therefore, **all files within a single analysis must have unique base filenames** to avoid collisions.
  - No spaces or special characters in file names
  - Case-sensitive matching between metadata and actual files
- **Location**: Files must be accessible from the workflow execution environment. `path_to_files_directory` is optional when absolute paths are provided in `fileName`.
- **Integrity**: Files should be properly formatted and not corrupted
- **Checksums**: `fileSize` and `fileMd5sum` are auto-calculated by the workflow; supply them in `file_metadata.tsv` only if you want the workflow to verify them against the actual files

### Metadata Files
- **Format**: Tab-separated values (TSV) with UTF-8 encoding
- **Required columns**: Each metadata file must include all required columns as specified in schemas below
- **Consistency**: Analysis IDs must be consistent across all metadata files
- **Relationships**: Foreign key relationships must be maintained between entities

### Metadata File Schema Examples

#### `file_metadata.tsv` (Required)

| Column                | Description                        | Example                | Required |
|-----------------------|------------------------------------|------------------------|----------|
| submitter_analysis_id | Unique analysis identifier         | analysis_001           | Yes      |
| fileName              | File name, subdirectory-relative path, or absolute path of the data file. | sample001.bam or data/wgs/sample001.bam or /absolute/path/sample001.bam | Yes      |
| fileSize              | Size in bytes (auto-calculated if omitted) | 123456789       | No       |
| fileMd5sum            | MD5 checksum of the file (auto-calculated if omitted) | 1a2b3c4d5e6f... | No       |
| fileType              | Type of file (CRAM/BAM/VCF/BCF)    | BAM                    | Yes      |
| fileAccess            | File access level                  | controlled             | Yes      |
| dataType              | Type of data content               | Aligned Reads          | Yes      |

> **Note on `fileSize` and `fileMd5sum`**: These values are **automatically calculated** from the actual data files during payload generation if not provided. If supplied, they are verified against the calculated values and an error is raised on mismatch. Omitting them removes the need to pre-compute checksums manually.

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
| variant_calling_strategy   | Strategy used for variant calling  | Tumour normal                    | No       |
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

### Pre-submission Validation Checklist

**📁 File Validation:**
- [ ] **File Existence**: All files listed in `file_metadata.tsv` are accessible (via `path_to_files_directory`, subdirectory path, or absolute path)
- [ ] **Index Files**: Each data file has its corresponding index (.crai, .bai, .tbi, .csi) alongside
- [ ] **File Integrity**: If `fileMd5sum` / `fileSize` are provided in `file_metadata.tsv`, they will be verified; otherwise the workflow calculates them automatically
- [ ] **File Access**: All files are readable by the workflow execution environment

**📊 Metadata Validation:**
- [ ] **Format Check**: TSV files open correctly in spreadsheet software without encoding issues
- [ ] **Required Columns**: All required columns are present with correct names
- [ ] **Data Consistency**: Analysis IDs are consistent across all metadata files
- [ ] **Controlled Vocabularies**: Values use permissible values

**🔧 Technical Validation:**
- [ ] **No Special Characters**: File names contain only alphanumeric characters, hyphens, and underscores
- [ ] **Path Accessibility**: All file paths are accessible from the workflow execution environment
- [ ] **Case Sensitivity**: File names in metadata exactly match actual file names

> 💡 **Quick Validation**: The workflow automatically performs these validations, but checking beforehand prevents submission failures.

### Generating MD5 Checksums
Pre-computing checksums is **no longer required**. The workflow automatically calculates `fileSize` and `fileMd5sum` for each file during the payload generation step. If you do provide these values in `file_metadata.tsv`, they will be verified against the calculated values; a mismatch will cause the analysis to fail.

If you still want to pre-compute checksums (e.g., for an audit trail before submission):
```bash
# For single file
md5sum sample001.cram

# For all files in directory
find /path/to/data -name "*.cram" -exec md5sum {} \; > checksums.txt
```

### Testing Your Input Setup
- Start with 1-2 small files to validate your environment and input setup
- See [Testing Guide](testing.md) for comprehensive testing guide using provided test dataset

### Recommended Directory Structure

The workflow supports several file organization patterns. Choose the one that best fits your setup.

#### Option A — Flat directory
All data files in a single directory, `fileName` contains the file name only:
```
study/
├── metadata/
│   ├── file_metadata.tsv          # fileName column: "sample1.cram"
│   └── analysis_metadata.tsv
└── data/                          # --path_to_files_directory
    ├── sample1.cram
    ├── sample1.cram.crai
    ├── sample2.cram
    └── sample2.cram.crai
```

#### Option B — Subdirectory pathing
Data files organised in subdirectories; `fileName` contains the relative sub-path from `path_to_files_directory`:
```
study/
├── metadata/
│   ├── file_metadata.tsv          # fileName column: "wgs/sample1.cram"
│   └── analysis_metadata.tsv
└── data/                          # --path_to_files_directory
    └── wgs/
        ├── sample1.cram
        ├── sample1.cram.crai
        ├── sample2.cram
        └── sample2.cram.crai
```

#### Option C — Absolute paths (no `path_to_files_directory`)
`fileName` contains an absolute path to each file; `--path_to_files_directory` can be omitted:
```
study/
└── metadata/
    ├── file_metadata.tsv          # fileName column: "/mnt/storage/data/sample1.cram"
    └── analysis_metadata.tsv
```

**Benefits of flexible pathing:**
- **Flat directory**: Simplest setup when all files reside in one place
- **Subdirectory pathing**: Organise large batches hierarchically without changing `path_to_files_directory`
- **Absolute pathing**: Ideal when files are spread across different storage systems or mount points

---

For further details, see:
- [Introduction Guide](introduction.md) - Workflow overview and architecture
- [Usage Guide](usage.md) - Step-by-step execution instructions
- [Testing Guide](testing.md) - Testing instructions with sample data
- [Receipt Guide](receipt.md) - Understanding submission results
- [Output Documentation](output.md) - Complete output file descriptions
- [Troubleshooting Guide (TBD)](troubleshooting.md) - Common issues and solutions