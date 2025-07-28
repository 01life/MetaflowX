process SUMMARYRESULT {
    
    label 'process_low'

    input:
    path(dastool_eval)

    output:
    path("${params.pipeline_prefix}_*")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    echo -e "sampleID\\tbinner\\tbin\\tbin_set\\tunique_SCGs\\tredundant_SCGs\\tSCG_set\\tsize\\tcontigs\\tN50\\tbin_score\\tSCG_completeness\\tSCG_redundancy" > header.xls
    
    cat header.xls ${dastool_eval} > all_eval.xls

    summaryDASToolsPlot.R all_eval.xls ${params.pipeline_prefix}

    """
    stub:
    """
    touch ${params.pipeline_prefix}_1.pdf
    touch ${params.pipeline_prefix}_1.png
    touch ${params.pipeline_prefix}_1.xls
    """
}



