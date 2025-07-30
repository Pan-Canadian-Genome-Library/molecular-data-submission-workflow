// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//https://github.com/Pan-Canadian-Genome-Library/Roadmap/issues/59
// DUMMY IMPLEMENTATION FOR TESTING - Replace with actual implementation

workflow CLINICAL_SERVICE_DATA_SUBMISSION {

    take:
    clinical_upload
    
    main:
    
    ch_versions = Channel.empty()
    
    // Create dummy versions
    ch_versions = ch_versions.mix(Channel.value("clinical_service: v1.0"))
    
    // Generate dummy successful registrations
    successful_registeration = clinical_upload
        .filter { meta, _csv_files -> meta.id != "analysis_002" } // Simulate that analysis_002 clinical data failed
        .map { meta, csv_files -> 
            // Add success metadata
            def success_meta = meta.clone()
            success_meta.registration_status = "success"
            success_meta.registered_id = "REG_${meta.id}"
            [success_meta, csv_files]
        }
    
    // Generate dummy unsuccessful registrations  
    unsuccessful_registeration = clinical_upload
        .filter { meta, _csv_files -> meta.id == "analysis_002" } // Simulate that analysis_002 clinical data failed
        .map { meta, csv_files ->
            // Add failure metadata
            def failure_meta = meta.clone() 
            failure_meta.registration_status = "failed"
            failure_meta.error_reason = "validation_failed"
            [failure_meta, csv_files]
        }
    
    emit:
    successful_registeration   // channel: [ val(meta), [csv] ] multiple CSVs per entity, each CSV contains unique identifier column
    unsuccessful_registeration // channel: [ val(meta), [csv] ] multiple CSVs per entity, each CSV contains unique identifier column  
    versions = ch_versions     // channel: [ versions.yml ]
}

