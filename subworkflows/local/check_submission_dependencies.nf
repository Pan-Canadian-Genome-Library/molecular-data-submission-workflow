// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/58

workflow CHECK_SUBMISSION_DEPENDENCIES {

    take:
        file_metadata // Spreadsheet
        analysis_metadata // Spreadsheet
        workflow_metadata // Spreadsheet
        read_group_metadata // Spreadsheet
        experiment_metadata // Spreadsheet
        specimen_metadata // Spreadsheet
        sample_metadata // Spreadsheet

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow
    Channel.value(1).subscribe{println "CHECK_SUBMISSION_DEPENDENCIES helloA"}

    clinical_upload = [[{},null]]
    entity_mapping = [[],{}]
    unsuccessful_dependency = [[],{}]
    molecular_files_to_upload = [[],{}]

    emit:
    // TODO nf-core: edit emitted channels
    clinical_upload  // channel: [ val(meta), [csv] ] multiple CSVs per entity
    molecular_files_to_upload // channel [ val(meta) [files] ] per analysis
    entity_mapping         // channel: [ val(meta), [ csv ] ] relational mapping
    unsuccessful_dependency // [val(meta),[csv]]
    versions = ch_versions                     // channel: [ versions.yml ]
}

