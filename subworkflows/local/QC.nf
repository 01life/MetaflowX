//
// QC
//


include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'
include { TRIM } from '../../modules/local/qc/trim'
include { FASTP } from '../../modules/local/qc/fastp'
include { FASTQSTAT } from '../../modules/local/qc/fastq_stat'
include { MERGEQC } from '../../modules/local/qc/merge_qc'
include { PIPELINEERROR } from '../../modules/local/common/pipeline_error'
include { PIPELINEEXIT } from '../../modules/local/common/pipeline_exit'

workflow QC {
    take:
    sample_number
    raw_reads

    main:

    /*
    * Verify the essential parameters for running this module
    */
    
    if (!(params.qc_tool in ["fastp", "trimmomatic"])) { exit 1, "The parameter qc_tool is invalid, supported values are:\n * fastp \n * trimmomatic" }

    qc_essential_db = [params.host_db, params.adapters]
    if(!checkEssentialParams(qc_essential_db)) { exit 1, "The required parameters to execute the QC module are:\n --host_db\n --adapters" }

    ch_host_db = params2Channel(params.host_db)
    ch_adapters = params2Channel(params.adapters)
    ch_phix_db = params2Channel(params.phix_db)

    fastp_json = Channel.empty()

    if (params.qc_tool == "fastp") {
        FASTP(raw_reads, ch_host_db, ch_adapters,ch_phix_db)
        clean_reads = FASTP.out.clean_reads
        stat_info = FASTP.out.stat_info
        fastp_json = FASTP.out.fastp_json.collect()
    } 
    
    if (params.qc_tool == "trimmomatic") {
        TRIM(raw_reads, ch_host_db, ch_adapters,ch_phix_db)
        clean_reads = TRIM.out.clean_reads
        stat_info = TRIM.out.stat_info
    }

    data = stat_info.join(clean_reads)
    
    FASTQSTAT(data)
    ch_reads_stat = FASTQSTAT.out.reads_stat.collect()
    finish_number = ch_reads_stat.flatten().filter { file -> !file.isEmpty() }.count()

    // QC completed successfully
    MERGEQC(sample_number, finish_number, ch_reads_stat, FASTQSTAT.out.clean_reads.collect())

    //generate an error log and terminate the pipeline if QC error occurs
    PIPELINEERROR("QC", sample_number, finish_number)
    PIPELINEEXIT(PIPELINEERROR.out.log)

    qc_report = MERGEQC.out.qc_report_list.mix(PIPELINEERROR.out.log).collect()

    emit:
    clean_reads // channel: [ val(id), [ reads1, reads2 ] ]
    qc_report
    fastp_json
    
}

