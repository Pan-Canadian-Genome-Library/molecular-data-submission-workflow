// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/58
include {CHECK_CLINICAL} from '../../modules/local/checkclinical'

workflow CHECK_SUBMISSION_DEPENDENCIES {

    take:
        study_id // Study ID
        token //Token for checking
        file_metadata // Spreadsheet
        analysis_metadata // Spreadsheet
        workflow_metadata // Spreadsheet
        read_group_metadata // Spreadsheet
        experiment_metadata // Spreadsheet
        specimen_metadata // Spreadsheet
        sample_metadata // Spreadsheet
        path_to_files_directory // file path

    main:


    ch_versions = Channel.empty()
    // Check clinical files
    CHECK_CLINICAL(
        study_id, // tuple(val : Study ID)
        token, // tuple(val : Token for checking)
        file_metadata, // tuple(val : Spreadsheet)
        analysis_metadata, // tuple(val : Spreadsheet)
        workflow_metadata, // tuple(val : Spreadsheet)
        read_group_metadata, // tuple(val : Spreadsheet)
        experiment_metadata, // tuple(val : Spreadsheet)
        specimen_metadata, // tuple(val : Spreadsheet)
        sample_metadata, // tuple(val : Spreadsheet)
        path_to_files_directory   //tuple(Path : directory)
    )
    // return analyses as channels
    CHECK_CLINICAL.out.status.flatten().collate(1).map{
        it ->
        [
            meta : [
                id : "${it[0].parent.toString().split('/')[-1]}",
                study : "${it[0].parent.toString().split('/')[-2]}",
            ],
            analysis : [
                analysis : file("${file(it[0]).parent}/*_analysis.tsv")[0],
                workflow : file("${file(it[0]).parent}/*_workflow.tsv")[0],
                files : file("${file(it[0]).parent}/*_files.tsv")[0]
            ],
            clinical : [
                specimen : file("${file(it[0]).parent}/*_specimen.tsv")[0],
                sample : file("${file(it[0]).parent}/*_sample.tsv")[0],
                experiment : file("${file(it[0]).parent}/*_experiment.tsv")[0],
                read_group : file("${file(it[0]).parent}/*_read_group.tsv")[0]
            ],
            status : it[0].toString().endsWith("SUCCESS"),
            status_file : file(it[0])

        ] 
    }.map{
        it -> {
            [
                meta : [
                    id : it.meta.id,
                    study : it.meta.study,
                    type : it.analysis.analysis.splitCsv(sep:'\t',header:true).analysisType[0]
                ],
                analysis : it.analysis,
                clinical : it.clinical,
                files : it.status ? it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { file("${path_to_files_directory[0]}/$it") } : [],
                status : it.status,
                status_file : it.status_file
            ]
        }
    }.set{analysis_channels}


    emit:
    // TODO nf-core: edit emitted channels
    analysis_channels // channel : [ 
    //    meta: [val(meta.id),val(meta.study),val(meta.type)], 
    //    analysis:[analysis.workflow,analysis.analysis,analysis.files], 
    //    clinical:[clinical.specimen,clinical.sample,clinical.experiment,clinical.read_group]
    //    files: [...]
    //    status : true/false
    //    status_file : id.SUCESS/id.FAILURE
    //]
    versions = ch_versions                     // channel: [ versions.yml ]
}

