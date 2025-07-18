// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/63
workflow METADATA_PAYLOAD_GENERATION {

    take:
    file_metadata // Spreadsheet
    analysis_metadata // Spreadsheet
    workflow_metadata // Spreadsheet
    
    main:
    channel.value(1).subscribe{println "METADATA_PAYLOAD_GENERATION helloC"}
    ch_versions = Channel.empty()
    files_to_upload = Channel.of([[],"A",[]],[[],"B",[]],[[],"C",[]])
    emit:
    files_to_upload  // channel: [ val(meta), JSON payload, [genomic files] ] per payload 

    versions = ch_versions                     // channel: [ versions.yml ]
}

