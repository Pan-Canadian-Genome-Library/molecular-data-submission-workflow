// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/59
include {CLINICAL_SUBMISSION} from '../../modules/local/clinical_submission/main.nf'

workflow CLINICAL_SERVICE_DATA_SUBMISSION {

    take:
    analysis_channels
    
    main:

    ch_versions = Channel.empty()

    CLINICAL_SUBMISSION(
        analysis_channels.map{
            it ->
            [
                meta : it.meta,
                analysis : it.analysis,
                clinical : it.clinical ,
                files : it.files
            ]
        }
    )

    registering_files = CLINICAL_SUBMISSION.out.analysis_channels
        .join(CLINICAL_SUBMISSION.out.status_file, by: 0)
        .map{
            meta, analysis, clinical, files, status_file ->

            def status_content = status_file.text
            def status_value = status_content.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()

            updated_meta.status = (status_value == 'pass' && meta.status == 'pass') ? 'pass2' : 'failed'
            [meta:updated_meta, analysis:analysis, clinical:clinical, files:files]
        }

    successful_registeration  = registering_files.filter{it -> it.meta.status=='pass2'}
    unsuccessful_registeration = registering_files.filter{it -> it.meta.status=='failed'}
    ch_versions = ch_versions.mix(CLINICAL_SUBMISSION.out.versions)

    emit:
    // TODO nf-core: edit emitted channels
    registering_files // [tuple(meta), tuple(analysis), tuple(clinical), tuple(files) ]
    successful_registeration // [tuple(meta), tuple(analysis), tuple(clinical), tuple(files) ] where meta.status=='pass'
    unsuccessful_registeration // [tuple(meta), tuple(analysis), tuple(clinical), tuple(files) ] where meta.status=='failed'
    versions = ch_versions                     // channel: [ versions.yml ]
    //    meta: [val(meta.id),val(meta.study),val(meta.type),meta.status], 
    //    analysis:[analysis.workflow,analysis.analysis,analysis.files], 
    //    clinical:[clinical.specimen,clinical.sample,clinical.experiment,clinical.read_group]
    //    files: [...]
}

