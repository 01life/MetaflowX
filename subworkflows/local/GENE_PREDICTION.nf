//
// gene prediction
//

include { PRODIGAL } from '../../modules/local/geneset/prodigal'

include { PIPELINEERROR } from '../../modules/local/common/pipeline_error'
include { PIPELINEEXIT } from '../../modules/local/common/pipeline_exit'

workflow GENEPREDICTION {
    take:
    sample_number
    contigs           // channel: [ val(id), path(contigs) ]

    main:

    prodigal_log = Channel.empty()

    if(params.prodigal_output){
        fpredict = file(params.prodigal_output, checkIfExists: true)
        // parse csv 
        predict_data = Channel.from (fpredict)
            .splitCsv (header:true, sep:',')
            .map { row ->
                // check the expected columns
                def expectedColumns = ['id', 'pep', 'cds'] 
                def missingColumns = expectedColumns.findAll { !row.containsKey(it) }
                if (missingColumns.size() > 0) {
                    exit 1, "Configuration error: Missing required columns in '--prodigal_output' parameter, verify: ${missingColumns.join(', ')}"
                }
                // check empty entries in expectedColumns
                for (column in expectedColumns) {  
                    if (row[column] == null || row[column].trim() == '') {  
                        exit(1, "Empty value found in column: ${column}")  
                    }  
                }  
                // get information
                return [row.id, row.pep, row.cds]  
            }
        prodigal_faa = predict_data.map{ it -> [it[0], it[1]]}
        prodigal_cds = predict_data.map{ it -> [it[0], it[2]]}
        // prodigal_noid_faa = prodigal_faa.map{ it -> it[1] }.collect()
        // prodigal_noid_cds = prodigal_cds.map{ it -> it[1] }.collect()
    // not run geneset prediction
    }else{
        PRODIGAL(contigs)
        prodigal_faa = PRODIGAL.out.faa     // channel: [ val(id),path(faa) ]
        prodigal_cds = PRODIGAL.out.cds
        prodigal_noid_faa = PRODIGAL.out.noid_faa.collect()
        // prodigal_noid_cds = PRODIGAL.out.noid_cds.collect()
        
        finish_number_prodigal = prodigal_noid_faa.flatten().filter { file -> !file.isEmpty() }.count()

        // generate an error log and terminate the pipeline if Prodigal error occurs
        PIPELINEERROR("Prodigal", sample_number, finish_number_prodigal)
        prodigal_log = PIPELINEERROR.out.log
        PIPELINEEXIT(prodigal_log)
    }

    emit:
    prodigal_faa
    prodigal_cds
    // prodigal_noid_faa
    // prodigal_noid_cds
    prodigal_log
}

