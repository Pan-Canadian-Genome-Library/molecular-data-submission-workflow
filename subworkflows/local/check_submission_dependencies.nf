// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/58
<<<<<<< Updated upstream
=======
include {CHECK_DEPENDENCIES} from '../../modules/local/check_dependencies'
include {ANALYSIS_SPLIT} from '../../modules/local/analysis_split'
include {VALIDATE_CLINICAL} from '../../modules/local/validate_clinical'
include {CLINICAL_SUBMISSION} from '../../modules/local/submit_clinical'
>>>>>>> Stashed changes

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
        // Check analysis and clinical files to ensure dependencies are met. Hard stop 
        CHECK_DEPENDENCIES(
            study_id, // tuple(val : Study ID)
            file_metadata, // tuple(val : Spreadsheet)
            analysis_metadata, // tuple(val : Spreadsheet)
            workflow_metadata, // tuple(val : Spreadsheet)
            read_group_metadata, // tuple(val : Spreadsheet)
            experiment_metadata, // tuple(val : Spreadsheet)
            specimen_metadata, // tuple(val : Spreadsheet)
            sample_metadata, // tuple(val : Spreadsheet)
            path_to_files_directory   //tuple(Path : directory)
        )
        ch_versions = ch_versions.mix(CHECK_DEPENDENCIES.out.versions)

<<<<<<< Updated upstream
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
=======
        //Split files (analysis and clinical) per analysis. Soft stop.
        ANALYSIS_SPLIT(
            study_id, // tuple(val : Study ID)
            file_metadata, // tuple(val : Spreadsheet)
            analysis_metadata, // tuple(val : Spreadsheet)
            workflow_metadata, // tuple(val : Spreadsheet)
            read_group_metadata, // tuple(val : Spreadsheet)
            experiment_metadata, // tuple(val : Spreadsheet)
            specimen_metadata, // tuple(val : Spreadsheet)
            sample_metadata, // tuple(val : Spreadsheet)
            CHECK_DEPENDENCIES.out.relational_mapping
        )
        ch_versions = ch_versions.mix(ANALYSIS_SPLIT.out.versions)
        //Using status file, reform channel to include meta, analysis, clinical, files, CHECK_DEPENDCIES output relational mapping, CHECK_DEPENDCIES output analysis_type, data directory path
        ANALYSIS_SPLIT.out.status.flatten()
        .combine(CHECK_DEPENDENCIES.out.relational_mapping)
        .combine(CHECK_DEPENDENCIES.out.analysis_types)
        .combine(path_to_files_directory)
        .map{
            it,relational_mapping,analysis_types,data_directory ->
            [
                meta : [
                    id : "${it.parent.toString().split('/')[-1]}",
                    study : "${it.parent.toString().split('/')[-2]}",
                    status : it.text.contains('status: "FAILED"') ? 'failed' : 'pass'
                ],
                analysis : [
                    analysis : file("${file(it).parent}/*analysis.tsv")[0],
                    workflow : file("${file(it).parent}/*workflow.tsv")[0],
                    files : file("${file(it).parent}/*files.tsv")[0]
                ],
                clinical : [
                    specimen : file("${file(it).parent}/*specimen.tsv")[0],
                    sample : file("${file(it).parent}/*sample.tsv")[0],
                    experiment : file("${file(it).parent}/*experiment.tsv")[0],
                    read_group : file("${file(it).parent}/*read_group.tsv")[0]
                ],
                status_file : file(it),
                relational_mapping: relational_mapping,
                analysis_types: analysis_types,
                data_directory : data_directory
            ] 
        }.map{ //update with file directory paths after, easier this way
            it -> {
                [
                    meta : [
                        id : it.meta.id,
                        study : it.meta.study,
                        type : workflow.stubRun ? "sequencingExperiment" : it.analysis.analysis.splitCsv(sep:'\t',header:true).analysisType[0],
                        status : it.meta.status
                    ],
                    analysis : it.analysis,
                    clinical : it.clinical,
                    files : workflow.stubRun ? [] : it.meta.status=='pass' ? it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { file("${path_to_files_directory[0]}/$it") } : [],
                    status_file : it.status_file,
                    relational_mapping: it.relational_mapping,
                    analysis_types: it.analysis_types,
                    data_directory: it.data_directory
                ]
            }
        }.set{analysis_channels}

        //files : it.meta.status=='pass' ? workflow.stubRun ? [] : it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { file("${path_to_files_directory[0]}/$it") } : [],
        //Validate clinical and analysis 
        VALIDATE_CLINICAL(
             analysis_channels
        )
        ch_versions = ch_versions.mix(VALIDATE_CLINICAL.out.versions)

        //Update channel status based on VALIDATE_CLINICAL outcome
        VALIDATE_CLINICAL.out.status.map{
            meta, analysis, clinical , files, status_file ->
            def status_value = status_file.text.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()

            updated_meta.status = (status_value == 'pass' && meta.status == 'pass') ? 'pass' : 'failed'
            [
                meta : updated_meta,
                analysis : analysis,
                clinical : clinical,
                files : files
            ]
        }.set{validated_analysis_channels}

        //Submit data to clinical.submission
        CLINICAL_SUBMISSION(validated_analysis_channels)

        //Update channel status based on CLINICAL_SUBMISSION outcome
         CLINICAL_SUBMISSION.out.status.map{
            meta, analysis, clinical , files, status_file ->
            def status_value = status_file.text.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()

            updated_meta.status = (status_value == 'pass' && meta.status == 'pass') ? 'pass' : 'failed'
            [
                meta : updated_meta,
                analysis : analysis,
                clinical : clinical,
                files : files
            ]
        }.set{submitted_analysis_channels}

    emit:
        versions = ch_versions
        analysis_channels = analysis_channels//format checked
        validated_analysis_channels = validated_analysis_channels //format checked and validated
        submitted_analysis_channels = submitted_analysis_channels//format checked, validated and submitted to clinical service
>>>>>>> Stashed changes
}

