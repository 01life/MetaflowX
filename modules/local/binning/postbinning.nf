process POSTBINNING {

    label 'process_low'

    input:
    path(final_bins)
    path(qs_quality_report)

    output:
    path("050.HQRawBin/HQBin/*.fa"),emit:"final_bins"
    path("050.HQRawBin/${params.pipeline_prefix}_all_Original_Bins_all_level_quality.xls"),emit:"qs_quality_report"
    path("binQS_report"),emit:"rawbinQS"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    cat ${qs_quality_report} > tmp.txt

    mkdir -p 050.HQRawBin/HQBin

    mv ${final_bins} 050.HQRawBin/HQBin
    
    cat tmp.txt | awk '!a[\$0]++' | uniq > 050.HQRawBin/${params.pipeline_prefix}_all_Original_Bins_all_level_quality.xls

    mkdir binQS_report 
    cp 050.HQRawBin/${params.pipeline_prefix}_all_Original_Bins_all_level_quality.xls binQS_report/rawbinQS.xls

    """

}