# Molecular Data submission workflow tests

## Test Cases

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
- checks subworkflow `heck_submission_dependencies`
- bad path demonstrates results when multiple analyses failed due to various conditions
- ensures subworkflow handles the analysis correctly