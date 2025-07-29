// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/58
// DUMMY IMPLEMENTATION FOR TESTING - Replace with actual implementation

include { PREPROCESS_SUBMISSION } from '../../../modules/local/submission/preprocess/main'

workflow CHECK_SUBMISSION_DEPENDENCIES {

    take:
        file_metadata // Spreadsheet
        analysis_metadata // Spreadsheet
        workflow_metadata // Spreadsheet
        read_group_metadata // Spreadsheet
        experiment_metadata // Spreadsheet
        specimen_metadata // Spreadsheet
        sample_metadata // Spreadsheet
        path_to_files_directory // file path to directory containing target files

    main:

    ch_versions = Channel.empty()
    
    // Process input metadata files to extract analysis information and create input channels
    // Read analysis metadata to get analysis IDs and types
    ch_input_data = Channel.fromPath(analysis_metadata)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [
                id: row.submitter_analysis_id,
                study: row.studyId,
                type: row.analysisType,
                status: "pass"
            ]
            [meta, file(file_metadata), file(analysis_metadata), file(workflow_metadata), file(read_group_metadata), file(experiment_metadata), file(specimen_metadata), file(sample_metadata), file(path_to_files_directory)]
        }
        .view { "ch_input_data: $it" }  // Add this line to view channel content


    // Use PREPROCESS_SUBMISSION module to process and reorganize the input files
    PREPROCESS_SUBMISSION(ch_input_data)
    ch_versions = ch_versions.mix(PREPROCESS_SUBMISSION.out.versions)

    // Process the generated data to create the expected output channels
    
    // Molecular files to upload - combine all metadata files with data files
    molecular_files_to_upload = PREPROCESS_SUBMISSION.out.file_meta
        .join(PREPROCESS_SUBMISSION.out.analysis_meta, by: 0)
        .join(PREPROCESS_SUBMISSION.out.workflow_meta, by: 0)
        .join(PREPROCESS_SUBMISSION.out.data_files, by: 0)
        .map { meta, file_meta, analysis_meta, workflow_meta, data_files_txt ->
            // Read the data files list from the text file
            def data_files_list = []
            if (data_files_txt.text.trim()) {
                data_files_list = data_files_txt.text.split('\n').findAll { it.trim() }.collect { file(it.trim()) }
            }
            [meta, file_meta, analysis_meta, workflow_meta, data_files_list]
        }

    // Biospecimen entity data - combine all entity files
    biospecimen_entity = PREPROCESS_SUBMISSION.out.specimen
        .join(PREPROCESS_SUBMISSION.out.sample, by: 0)
        .join(PREPROCESS_SUBMISSION.out.experiment, by: 0)
        .join(PREPROCESS_SUBMISSION.out.read_group, by: 0)
        .map { meta, specimen, sample_file, experiment, read_group ->
            [meta, [specimen, sample_file, experiment, read_group]]
        }

    // Status files
    status = PREPROCESS_SUBMISSION.out.status

    emit:
    molecular_files_to_upload     // channel: [ val(meta), path(file_meta.tsv), path(analysis_meta.tsv), path(workflow_meta.tsv), [data_files] ] per analysis
    biospecimen_entity               // channel: [ val(meta), [ specimen, sample, experiment, read_group ] ] relational mapping
    status                       // channel: [ val(meta), path(*_status.yml) ]
    versions = ch_versions       // channel: [ versions.yml ]
}

