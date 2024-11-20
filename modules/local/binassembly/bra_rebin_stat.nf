
process REBINSTAT {

    tag "$type"

    label 'process_low'

    input:
    val(type)
    path(reas_bins)

    output:
    tuple val("reassembly_bins"),path("Reassembly_Bins"),emit:"reassembly_bins"
    path("${params.pipeline_prefix}_${type}_bins_info.xls"),emit:"contig_info"
    // path("REBINSTAT_report") ,emit:"report"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    mkdir Reassembly_Bins
    mv ${reas_bins} Reassembly_Bins/

    ls Reassembly_Bins/*fa* |awk -F "/" '{print\$NF"\\t"\$0}' |sed 's/_contigs.fa\\t/\\t/g' > bin.fa.list

    bin_stat.py bin.fa.list ${params.pipeline_prefix}_${type}_bins_info.xls

    #mkdir REBINSTAT_report
    #cp ${params.pipeline_prefix}_${type}_bins_info.xls REBINSTAT_report/reasbinInfo.xls

    """

}
