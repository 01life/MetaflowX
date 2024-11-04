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

    ls ${dastool_eval} | while read a ;do tail -n +2 \$a  ;done > all_eval.xls

    sed -i '1i sampleID\\tbinner\\tbin\\tbin_set\\tunique_SCGs\\tredundant_SCGs\\tSCG_set\\tsize\\tcontigs\\tN50\\tbin_score\\tSCG_completeness\\tSCG_redundancy' all_eval.xls

    summaryDASToolsPlot.R all_eval.xls ${params.pipeline_prefix}

    """
    
}



