// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


workflow SUBMISSION_RECEIPT {

    take:
       unsuccessful_dependency // [val(meta),[csv]]
       unsuccessful_registeration // channel: [ val(meta), [csv] ] multiple CSVs per entity
       unsuccessful_validation  //[ val(meta), JSON payload, [genomic files] , metadata sanity boolean, sequencing integrity boolean , cross validation boolean]
       unsuccessful_upload
       successful_registeration // channel: [ val(meta), [csv] ] multiple CSVs per entity
       successful_validation  //[ val(meta), JSON payload, [genomic files] , metadata sanity boolean, sequencing integrity boolean , cross validation boolean]
       successful_upload
    
    main:

    ch_versions = Channel.empty()
    Channel.value(1).subscribe{println "SUBMISSION_RECEIPT helloE"}

    emit:
    // TODO nf-core: edit emitted channels

    versions = ch_versions                     // channel: [ versions.yml ]
}

