# Molecular Data submission workflow tests

## Test Cases

### Workflows

### Subworkflows
#### 1. Test Subworkflow CHECK_SUBMISSION_DEPENDENCIES
1. check check_submission_dependencies happy path
- > [!IMPORTANT]<Br>Before running tests, ensure that the correct API token is provided through `tests/nextflow.config`. Otherwise `happy path` tests will fail at earlier
- > Before testing, ensure only the entity present is participant. Otherwise the system will detect previously submtited non-participant entities and submission will be bypassed.<Br>
  - > To do so, the example study used is `EXAMPLE-CA9`. Query `specimen` data for `NEW_SPECIMEN_02` in study `EXAMPLE-CA9` using the `/data/category/{categoryId}/organization/{organization}`. Use the retrieve `specimen` `systemId` as the arguements for the DELETE endpoint `/submission/category/{categoryId}/data/{systemId}`.<Br>
  - > This should generate a `submissionId` to be used in the commit POST endpoint `/submission/category/{categoryId}/commit/{submissionId}`. Recheck data to ensure only donors exist.<Br>
  - > Note for `categoryId`: unless a category is named after the organization, defaults to using `prod_pcgl_schema`
- the following shoould be the results for bulk submission
  - `NEW_ANALYSIS_01` succeeds
  - `NEW_ANALYSIS_02` fails b/c of inconsistent records
  - `DONOR_01_EXPERIMENT_02` fails b/c of missing foreign key
2. check check_submission_dependencies stub
 
### Modules

#### 1. Test module CHECK_DEPENDENCIES
1. check dependencies happy path - success
- > [!IMPORTANT]<Br>Before running tests, ensure that the correct API token is provided through `tests/nextflow.config`. Otherwise `happy path` tests will fail at earlier
2. check dependencies bad path - fail
- fails b/c of wrong study, returning `404`
3. check dependencies happy path stub

#### 2. Test module ANALYSIS_SPLIT
1. analysis split happy path - success
2. analysis split no biospecimens - success
3. analysis split bad path - fail
- 1/4 analysis should fail b/c of conflicting analysis records : `NON_EXISTING_EXPERIMENT` fails b/c of conflicting analysis records
4. analysis split stub

#### 3. Test module VALIDATE_CLINICAL
1. Validate clinical happy path
- > [!IMPORTANT]<Br>Before running tests, ensure that the correct API token is provided through `tests/nextflow.config`. Otherwise `happy path` tests will fail at earlier
2. Validate clinical bad path - fail dependency
- fails b/c prior status is fail
3. Validate clinical bad path - fail submission
- > [!IMPORTANT]<Br>Before running tests, ensure that the correct API token is provided through `tests/nextflow.config`. Otherwise `happy path` tests will fail at earlier
- fails for the following reasons:
  - Field 'experiment_design' is not consistent for record NEW_EXPERIMENT_02 in entity experiment. Specified - 4 rounds of PCR vs Comitted - 6 rounds of PCR
  - Field 'submitter_sample_id' is not consistent for record NEW_EXPERIMENT_02 in entity experiment. Specified - NEW_SAMPLE_02 vs Comitted - ERROR_IN_SAMPLE
  - File bad_path/D0RE2.1_r1.fq.gz could not be found
4. Validate clinical stub

#### 4. Test module CLINICAL_SUBMISSION
1. clinical service data submission happy path
- > [!IMPORTANT]<Br>Before running tests, ensure that the correct API token is provided through `tests/nextflow.config`. Otherwise `happy path` tests will fail at earlier
- > Before testing, ensure only the entity present is participant. Otherwise the system will detect previously submtited non-participant entities and submission will be bypassed.<Br>
  - > To do so, the example study used is `EXAMPLE-CA9`. Query `specimen` data for `NEW_SPECIMEN_02` in study `EXAMPLE-CA9` using the `/data/category/{categoryId}/organization/{organization}`. Use the retrieve `specimen` `systemId` as the arguements for the DELETE endpoint `/submission/category/{categoryId}/data/{systemId}`.<Br>
  - > This should generate a `submissionId` to be used in the commit POST endpoint `/submission/category/{categoryId}/commit/{submissionId}`. Recheck data to ensure only donors exist.<Br>
  - > Note for `categoryId`: unless a category is named after the organization, defaults to using `prod_pcgl_schema`
2. clinical service data submission happy path - redundant
- > [!IMPORTANT]<Br>Before running tests, ensure that the correct API token is provided through `tests/nextflow.config`. Otherwise `happy path` tests will fail at earlier
- > Before testing, ensure only the entity present is participant. Otherwise the system will detect previously submtited non-participant entities and submission will be bypassed.<Br>
  - > To do so, the example study used is `EXAMPLE-CA9`. Query `specimen` data for `NONREDUNDANT_SAMPLE_02` in study `EXAMPLE-CA9` using the `/data/category/{categoryId}/organization/{organization}`. Use the retrieve `specimen` `systemId` as the arguements for the DELETE endpoint `/submission/category/{categoryId}/data/{systemId}`.<Br>
  - > This should generate a `submissionId` to be used in the commit POST endpoint `/submission/category/{categoryId}/commit/{submissionId}`. Recheck data to ensure only donors exist.<Br>
  - > Note for `categoryId`: unless a category is named after the organization, defaults to using `prod_pcgl_schema`
- two records are submitted `NONREDUNDANT_SAMPLE_02` and `REDUNDANT_SPECIMEN_02`. Only `NONREDUNDANT_SAMPLE_02` is accepted b/c `REDUNDANT_SPECIMEN_02` was previously submitted
3. clinical service data submission bad path - fail dependency
4. clinical service data submission bad path - fail submission
- > [!IMPORTANT]<Br>Before running tests, ensure that the correct API token is provided through `tests/nextflow.config`. Otherwise `happy path` tests will fail at earlier∂
- fails for the following reasons:
  - invalid foreignKey: participants do not exist
  - results in downstream invalid foreign keys
  - duplicate read group records
5. clinical service data submission stub