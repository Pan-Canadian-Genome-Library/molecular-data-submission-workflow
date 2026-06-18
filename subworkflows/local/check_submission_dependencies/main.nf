// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/58

include {CHECK_DEPENDENCIES} from '../../../modules/local/check_dependencies'
include {ANALYSIS_SPLIT} from '../../../modules/local/analysis_split'
include {VALIDATE_CLINICAL} from '../../../modules/local/validate_clinical'
include {CLINICAL_SUBMISSION} from '../../../modules/local/submit_clinical'
include {LECTERN_VALIDATE} from '../../../modules/local/lectern/validate'

workflow CHECK_SUBMISSION_DEPENDENCIES {

    take:
    study_id // tuple (val:study_id)
    file_metadata // Spreadsheet
    analysis_metadata // Spreadsheet
    workflow_metadata // Spreadsheet
    read_group_metadata // Spreadsheet
    experiment_metadata // Spreadsheet
    specimen_metadata // Spreadsheet
    sample_metadata // Spreadsheet
    path_to_files_directory // file path

    main:
    
        ch_versions = channel.empty()
        // Check analysis and clinical files to ensure dependencies are met. Hard stop 
        CHECK_DEPENDENCIES(
            study_id, // tuple (val:study_id)
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

        // Check if minimum requirements are met - graceful stop if not
        // Update status channel with CHECK_DEPENDENCIES output
        CHECK_DEPENDENCIES.out.status
            .subscribe { _study, status_file ->
                def status_content = status_file.text
                def is_failed = status_content.contains('status: "FAILED"')
                
                if (is_failed) {
                    // Extract error details using functional approach
                    def lines = status_content.split('\n')
                    def error_start_index = lines.findIndexOf { item -> item.trim().startsWith('error_details:') }
                    def error_details = ""
                    
                    if (error_start_index >= 0) {
                        // Get lines after error_details: that start with spaces
                        def error_lines = lines[(error_start_index+1)..-1]
                            .takeWhile { item -> item.startsWith('    ') || item.trim().isEmpty() }
                            .collect { item -> item.replaceFirst('    ', '') }
                        
                        error_details = error_lines.join('\n').trim()
                    }                  
                    // Print user-friendly error message
                    System.err.println """
╔══════════════════════════════════════════════════════════════════════════════╗
║                        🚨 WORKFLOW STOPPED                                   ║
║              Minimum requirements for data submission not met!               ║
╚══════════════════════════════════════════════════════════════════════════════╝
""".stripIndent()
System.err.println """
❌ Study: ${study_id}

🔍 Issues found:
--------------------------------------------------------------------------------""".stripIndent()
System.err.println """${error_details}""".stripIndent()
System.err.println """--------------------------------------------------------------------------------

📋 Common solutions:
   • Ensure that access token with correct submission scope is provided
   • Ensure that the study has been registered
   • Ensure all required metadata files are provided and correctly formatted
   • Check that file paths in metadata files are valid and files exist
   • Verify that all required columns are present in metadata files
   • Ensure study_id matches across all metadata files
""".stripIndent()
System.err.println """
💡 For detailed error information, check: ${status_file}
""".stripIndent()
System.err.println """
Please fix the above issues and re-run the workflow.
                    """.stripIndent()
                    
                    // Flush the error stream before exiting
                    System.err.flush()

                    // Add a small delay to ensure output is displayed
                    Thread.sleep(00)

                    // Gracefully exit the workflow
                    System.exit(1)
                }
            }

        // Continue with the rest of the workflow only if dependencies check passed
        CHECK_DEPENDENCIES.out.status
            .filter { _study, status_file -> 
                !status_file.text.contains('status: "FAILED"')
            }
            .subscribe {
                log.info "✅ Minimum requirements check passed. Proceeding with data submission..."
            }

        //Split files (analysis and clinical) per analysis. Soft stop.
        ANALYSIS_SPLIT(
            study_id,
            file_metadata, // tuple(val : Spreadsheet)
            analysis_metadata, // tuple(val : Spreadsheet)
            workflow_metadata, // tuple(val : Spreadsheet)
            read_group_metadata, // tuple(val : Spreadsheet)
            experiment_metadata, // tuple(val : Spreadsheet)
            specimen_metadata, // tuple(val : Spreadsheet)
            sample_metadata, // tuple(val : Spreadsheet)
            CHECK_DEPENDENCIES.out.relational_mapping,
            params.path_to_files_directory ? //Mount provided data directory path or absolute file paths
                file_metadata.splitCsv(sep:'\t',header: true).filter { !it.fileName.startsWith('s3') }.map{ row -> file("${params.path_to_files_directory}/${row.fileName}").parent }.unique().collect() :\
                file_metadata.splitCsv(sep:'\t',header: true).filter { !it.fileName.startsWith('s3') }.map{ row -> file(row.fileName).parent }.unique().collect()
        )
        ch_versions = ch_versions.mix(ANALYSIS_SPLIT.out.versions)

        // Downstream receipt is only provided if join with files is successful. This check identifies the event where all analyses fail from start, usually b/c files cannot be found
        ANALYSIS_SPLIT.out.status.flatten().map{
            it ->
            [
                it.text.contains('status: "FAILED"') ? "failed" : "pass",
                "${it}"
            ]
        }
        .collect(flat:false)
        .transpose()
        .collect(flat:false)//.subscribe{println "${it}"}
        .map{
            it,_files ->
            // Pass as long as one analysis is successful in splitting
            def is_failed = it.join(',').contains('pass') ? false : true

            if (is_failed==true) {
                System.err.flush()

                System.err.println """
╔══════════════════════════════════════════════════════════════════════════════╗
║                        🚨 WORKFLOW STOPPED                                   ║
║                       No viable analyses found!                              ║
║                  Associated files could not be found                         ║
╚══════════════════════════════════════════════════════════════════════════════╝
""".stripIndent()
                    
                // Flush the error stream before exiting
                //System.err.flush()

                // Add a small delay to ensure output is displayed
                Thread.sleep(00)

                // Gracefully exit the workflow
                System.exit(1)
            }
        }
 

        if (params.path_to_files_directory){
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
                        analysis : files("${file(it).parent}/*analysis.tsv")[0],
                        workflow : files("${file(it).parent}/*workflow.tsv")[0],
                        files : files("${file(it).parent}/*files.tsv")[0]
                    ],
                    clinical : [
                        specimen : files("${file(it).parent}/*specimen.tsv")[0],
                        sample : files("${file(it).parent}/*sample.tsv")[0],
                        experiment : files("${file(it).parent}/*experiment.tsv")[0],
                        read_group : files("${file(it).parent}/*read_group.tsv")[0]
                    ],
                    status_file : file(it),
                    relational_mapping: relational_mapping,
                    analysis_types: analysis_types,
                    data_directory : data_directory
                ] 
            }.map{ //update with file directory paths after, easier this way
                it ->
                    [
                        meta : [
                            id : it.meta.id,
                            study : it.meta.study,
                            type : workflow.stubRun ? "sequencingExperiment" : it.analysis.analysis.splitCsv(sep:'\t',header:true).analysisType[0],
                            status : it.meta.status
                        ],
                        analysis : it.analysis,
                        clinical : it.clinical,
                        //| Condition                                               | Result                               |
                        //| ------------------------------------------------------- | ------------------------------------ |
                        //| `workflow.stubRun == true`                              | `[]`                                 |
                        //| `meta.status != 'pass'`                                 | `[]`                                 |
                        //| `meta.status == 'pass'` **and** `data_directory exists` | files with `data_directory/filename` |
                        //| `meta.status == 'pass'` **and** no `data_directory`     | files using filename only            |
                        //files : workflow.stubRun ? [] : it.meta.status=='pass' ? it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { fname -> file("${it.data_directory}/${fname}") } : [],
                        files : workflow.stubRun ? [] : it.meta.status=='pass' ? it.data_directory ? it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { fname -> file("${it.data_directory}/${fname}") } : it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { fname -> file("${fname}") }  : [],
                        status_file : it.status_file,
                        relational_mapping: it.relational_mapping,
                        analysis_types: it.analysis_types,
                        data_directory: it.data_directory
                    ]
               
            }.set{analysis_channels}
        } else {
            ANALYSIS_SPLIT.out.status.flatten()
            .combine(CHECK_DEPENDENCIES.out.relational_mapping)
            .combine(CHECK_DEPENDENCIES.out.analysis_types)
            .map{
                it,relational_mapping,analysis_types ->
                [
                    meta : [
                        id : "${it.parent.toString().split('/')[-1]}",
                        study : "${it.parent.toString().split('/')[-2]}",
                        status : it.text.contains('status: "FAILED"') ? 'failed' : 'pass'
                    ],
                    analysis : [
                        analysis : files("${file(it).parent}/*analysis.tsv")[0],
                        workflow : files("${file(it).parent}/*workflow.tsv")[0],
                        files : files("${file(it).parent}/*files.tsv")[0]
                    ],
                    clinical : [
                        specimen : files("${file(it).parent}/*specimen.tsv")[0],
                        sample : files("${file(it).parent}/*sample.tsv")[0],
                        experiment : files("${file(it).parent}/*experiment.tsv")[0],
                        read_group : files("${file(it).parent}/*read_group.tsv")[0]
                    ],
                    status_file : file(it),
                    relational_mapping: relational_mapping,
                    analysis_types: analysis_types,
                    data_directory : []
                ] 
            }.map{ //update with file directory paths after, easier this way
                it ->
                    [
                        meta : [
                            id : it.meta.id,
                            study : it.meta.study,
                            type : workflow.stubRun ? "sequencingExperiment" : it.analysis.analysis.splitCsv(sep:'\t',header:true).analysisType[0],
                            status : it.meta.status
                        ],
                        analysis : it.analysis,
                        clinical : it.clinical,
                        //| Condition                                               | Result                               |
                        //| ------------------------------------------------------- | ------------------------------------ |
                        //| `workflow.stubRun == true`                              | `[]`                                 |
                        //| `meta.status != 'pass'`                                 | `[]`                                 |
                        //| `meta.status == 'pass'` **and** `data_directory exists` | files with `data_directory/filename` |
                        //| `meta.status == 'pass'` **and** no `data_directory`     | files using filename only            |
                        //files : workflow.stubRun ? [] : it.meta.status=='pass' ? it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { fname -> file("${it.data_directory}/${fname}") } : [],
                        files : workflow.stubRun ? [] : it.meta.status=='pass' ? it.data_directory ? it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { fname -> file("${it.data_directory}/${fname}") } : it.analysis.files.splitCsv(sep:'\t',header:true).fileName.collect { fname -> file("${fname}") }  : [],
                        status_file : it.status_file,
                        relational_mapping: it.relational_mapping,
                        analysis_types: it.analysis_types,
                        data_directory: it.data_directory
                    ]
                
            }.set{analysis_channels}  
        }


        //Validate clinical and analysis required entities co-exist
        VALIDATE_CLINICAL(
             analysis_channels
        )
        ch_versions = ch_versions.mix(VALIDATE_CLINICAL.out.versions)

        //Update channel status based on VALIDATE_CLINICAL outcome
        VALIDATE_CLINICAL.out.status.map{
            meta, analysis, clinical , files, status_file, relational_mapping, analysis_types, data_directory ->
            def status_value = status_file.text.contains('status: "FAILED"') ? 'failed' : 'pass'
            def updated_meta = meta.clone()

            updated_meta.status = (status_value == 'pass' && meta.status == 'pass') ? 'pass' : 'failed'
            [
                meta : updated_meta,
                analysis : analysis,
                clinical : clinical,
                files : files,
                status_file : status_file,
                relational_mapping: relational_mapping,
                analysis_types: analysis_types,
                data_directory: data_directory
            ]
        }.set{validated_analysis_channels}

        ////validated_analysis_channels.subscribe{println "validated_analysis_channels: $it"}
        if (params.skip_upload){
            //validated_analysis_channels.subscribe{println "${it}"}

            LECTERN_VALIDATE(validated_analysis_channels)
            ////LECTERN_VALIDATE.out.status.subscribe{println "LECTERN_VALIDATE: $it"}
            //Update channel status based on CLINICAL_SUBMISSION outcome
            LECTERN_VALIDATE.out.status.map{
                meta, analysis, clinical , files, status_file ->
                def status_value = status_file.text.contains('status: "FAILED"') ? 'failed' : 'pass'
                def updated_meta = meta.clone()

                updated_meta.status = (status_value == 'pass' && meta.status == 'pass') ? 'pass' : 'failed'
                [
                    meta : updated_meta,
                    analysis : analysis,
                    clinical :  clinical,
                    files : files,
                    status_file : status_file
                ]
            }.set{submitted_analysis_channels}
        } else {
            //Submit data to clinical.submission
            CLINICAL_SUBMISSION(validated_analysis_channels)

            //Update channel status based on CLINICAL_SUBMISSION outcome
            CLINICAL_SUBMISSION.out.status.map{
                //Briefly rename files to genomic files to avoid collision with files()
                meta, analysis, clinical , genomic_files, status_file ->
                def status_value = status_file.text.contains('status: "FAILED"') ? 'failed' : 'pass'
                def updated_meta = meta.clone()

                updated_meta.status = (status_value == 'pass' && meta.status == 'pass') ? 'pass' : 'failed'
                def parentDir = file(status_file).parent
                def specimen_files   = files("${parentDir}/*/retrieved/specimen.tsv")
                def sample_files     = files("${parentDir}/*/retrieved/sample.tsv")
                def experiment_files = files("${parentDir}/*/retrieved/experiment.tsv")
                def read_group_files  = files("${parentDir}/*/retrieved/read_group.tsv")

                def updated_clinical = [
                    specimen : specimen_files.size()==0 ? null :  specimen_files[0],
                    sample : sample_files.size()==0 ? null :  sample_files[0],
                    experiment : experiment_files==0 ? null :  experiment_files[0],
                    read_group : read_group_files.size()==0 ? null : read_group_files[0]
                ]
                [
                    meta : updated_meta,
                    analysis : analysis,
                    clinical :  updated_clinical,
                    files : genomic_files,
                    status_file : status_file
                ]
            }.set{submitted_analysis_channels}          
        }


        analysis_channels.map{
            it ->
            [
                meta : it.meta,
                status_file : it.status_file
            ]
        }.collect().combine(
            validated_analysis_channels.map{
            it ->
                [
                meta : it.meta,
                status_file : it.status_file
                ]

            }.collect()
        ).combine(
            submitted_analysis_channels.map{
            it ->
                [
                meta : it.meta,
                status_file : it.status_file
                ]

            }.collect()
        ).set{status_files}

    emit:
        versions = ch_versions
        analysis_channels = analysis_channels//format checked
        validated_analysis_channels = validated_analysis_channels //format checked and validated
        submitted_analysis_channels = submitted_analysis_channels//format checked, validated and submitted to clinical service (only if skip upload is false)
        status = status_files //format checked, validated and submitted to clinical service
}
