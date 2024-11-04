//
// Processing of visualization reports and result files.
//

include { params2Channel } from '../../modules/local/common/utils'
include { REPORT } from '../../modules/local/polish/report'
include { POSTPIPELINE } from '../../modules/local/polish/postpipeline'
include { MULTIQC as MULTIQCFASTP } from '../../modules/local/polish/multiQC'
include { MULTIQC as MULTIQCBOWTIE2 } from '../../modules/local/polish/multiQC'

ch_report_topic = params2Channel(params.report_topic)
ch_report_order = params2Channel(params.report_order)
ch_report_template = params2Channel(params.report_template)
ch_report_images = params2Channel(params.report_images)

workflow POLISH {
    take:
    report
    fastp_json
    bowtie2_log
    clean_reads
    raw_reads_flag

    main:
    
    MULTIQCFASTP("fastp", fastp_json)
    
    MULTIQCBOWTIE2("bowtie2", bowtie2_log)

    REPORT(report, ch_report_topic, ch_report_order, ch_report_template, ch_report_images)

    //Process clean reads only when raw reads are input.
    if( params.mode==1 || (params.mode==0 && !params.skip_qc) || (params.mode in [2, 3] && raw_reads_flag) ){
        POSTPIPELINE(REPORT.out.html, clean_reads)
    }
    

}

