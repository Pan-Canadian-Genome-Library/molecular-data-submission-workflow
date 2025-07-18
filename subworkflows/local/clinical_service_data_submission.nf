// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


workflow CLINICAL_SERVICE_DATA_SUBMISSION {

    take:
    clinical_upload
    
    main:
    Channel.value(1).subscribe{println "CLINICAL_SERVICE_DATA_SUBMISSION helloB"}
    ch_versions = Channel.empty()
    successful_registeration = [{},[]]
    unsuccessful_registeration = [{},[]]
    emit:
    // TODO nf-core: edit emitted channels
    successful_registeration // channel: [ val(meta), [csv] ] multiple CSVs per entity, each CSV contains unique identifer column
    unsuccessful_registeration // channel: [ val(meta), [csv] ] multiple CSVs per entity, each CSV contains unique identifer column
    versions = ch_versions                     // channel: [ versions.yml ]
}

