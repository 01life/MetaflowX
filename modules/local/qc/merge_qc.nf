
process MERGEQC {

    label 'process_single'

    input:
    val(sample_number)
    val(finish_number)
    path(reads_stat)
    path(clean_reads)

    output:
    path("${params.pipeline_prefix}_all_sample_reads_stat.xls")
    path("qc_report"),emit:"qc_report_list"
    path("clean.reads.path.list")

    when:
    sample_number == finish_number

    script:
    """
    
    cat ${reads_stat} > all.stat.txt
    head -n 1 all.stat.txt > p0
    grep -v "ID" all.stat.txt > p1

    cat p0 p1 > ${params.pipeline_prefix}_all_sample_reads_stat.xls

    rm -rf p0 p1 all.stat.txt

    mkdir qc_report
    cp ${params.pipeline_prefix}_all_sample_reads_stat.xls qc_report/readstat.xls

    cat ${clean_reads} > clean.reads.path.list

    """

}
