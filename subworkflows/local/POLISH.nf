//
// Processing of visualization reports and result files.
//

include { params2Channel } from '../../modules/local/common/utils'
include { REPORT } from '../../modules/local/polish/report'
include { POSTPIPELINE } from '../../modules/local/polish/postpipeline'
include { MULTIQC as MULTIQCFASTP } from '../../modules/local/polish/multiQC'
include { MULTIQC as MULTIQCBOWTIE2 } from '../../modules/local/polish/multiQC'
include { NOTIFICATION } from '../../modules/local/common/notification'

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

    report.flatten().branch { file ->
        log: file.name.contains('log')
        result: true
    }
    .set { report_input }

    all_log = report_input.log.collect()
    all_report = report_input.result.collect()

    REPORT(all_report, all_log.ifEmpty([]), ch_report_topic, ch_report_order, ch_report_template, ch_report_images)

    //Process clean reads only when raw reads are input.
    if( params.mode==1 || (params.mode==0 && !params.skip_qc) || (params.mode in [2, 3] && raw_reads_flag) ){
        POSTPIPELINE(clean_reads.combine(REPORT.out.html))
    }
    
    //WeChat robot group notification: workflow execution completed.
    if(params.webhookurl){
        content = "🤩 Awesome! ${params.pipeline_prefix} has generated the report and is wrapping things up. Almost there! 👏"
        NOTIFICATION("POLISH", content, REPORT.out.html)
    }

}

