# Molecular Data submission workflow tests

## Test Cases

### check_dependencies

> [!IMPORTANT]
> Before running tests, ensure that the correct API token is provided through `tests/nextflow.config`. Otherwise `happy path` tests will fail at earlier unintended points 

1. Check clinical happy path
- checks module `check_clinical`
- happy path demonstrates result when all provided clinical files are correct
2. Check clinical bad path
- checks module `check_clinical`
- bad path demonstrates results when multiple analyses failed due to various conditions
- ensures workflows gracefully handles failed checks
3. Check check_submission_dependencies happy path
- checks subworkflow `heck_submission_dependencies`
- happy path demonstrates result when all provided clinical files are correct
- ensures subworkflow handles the analysis correctly
4. Check check_submission_dependencies bad path
- checks subworkflow `check_submission_dependencies`
- bad path demonstrates results when multiple analyses failed due to various conditions
- ensures subworkflow handles the analysis correctly

### Submit clinical

> [!IMPORTANT]
> Before testing, ensure only the entity present is participant. Otherwise the system will detect previously submtited non-participant entities and block the test from successful completion.<Br>
> To do so, the example study used is `EXAMPLE-CA9`. Query `specimen` data for `EXAMPLE-CA` using the `/data/category/{categoryId}/organization/{organization}`. Use the retrieve `specimen` `systemId` as the arguements for the DELETE endpoint `/submission/category/{categoryId}/data/{systemId}`.<Br>
> This should generate a `submissionId` to be used in the commit POST endpoint `/submission/category/{categoryId}/commit/{submissionId}`. Recheck data to ensure only donors exist.<Br>
> Note for `categoryId`: unless a category is named after the organization, defaults to using `prod_pcgl_schema`
