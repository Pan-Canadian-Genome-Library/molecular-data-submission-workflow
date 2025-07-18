// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/60
//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/61

workflow DATA_VALIDATION {

    take:
    files_to_upload //[ val(meta), JSON payload, [genomic files] ]
    
    main:
    Channel.value(1).subscribe{println "DATA_VALIDATION helloD"}
    ch_versions = Channel.empty()
    successful_validation = Channel.of([[],"A",[],true,true,true],[[],"C",[],true,true,true])
    unsuccessful_validation = Channel.of([[],"B",[],false,false,false])

    emit:
    // TODO nf-core: edit emitted cha nnels
    successful_validation //[ val(meta), JSON payload, [genomic files] , metadata sanity boolean, sequencing integrity boolean , cross validation boolean]
    unsuccessful_validation //[ val(meta), JSON payload, [genomic files]  , metadata sanity boolean, sequencing integrity boolean , cross validation boolean]

    versions = ch_versions                     // channel: [ versions.yml ]
}

