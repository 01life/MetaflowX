process MERGEQUAST {

    label 'process_low'
    
    input:
    path(reports)


    output:
    path("${params.pipeline_prefix}_quast_report.txt"),emit:"contig_stat"
    // path("${params.pipeline_prefix}_quast_report.html"),emit:"html_report"
    path("quast_report"),emit:"quast_report"

    script:
    """

    awk 'FNR==1 && NR!=1 {next;}{print}'  ${reports} > ${params.pipeline_prefix}_quast_report.txt
    
    #merge_and_plot_quast.py --result_dir ./ --prefix  ${params.pipeline_prefix}

    mkdir quast_report
    cp ${params.pipeline_prefix}_quast_report.txt quast_report/QuastContigstat.txt
    
    """

    stub:
    """
    touch ${params.pipeline_prefix}_quast_report.txt
    touch ${params.pipeline_prefix}_quast_report.html
    mkdir quast_report
    """
}
