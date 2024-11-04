process  PICKPBODBO {
    
    tag "$id"
    
    label 'process_low'

    input:
    val(optimizeBins_method)
    tuple val(id), path(bin_floder),path(bin_qs)

    output:
    tuple val(id),path("${id}_Best_Bin"),path("${id}_Best_Bin_quality_report.tsv"),emit:"best"

    when:
    task.ext.when == null || task.ext.when

    script:
    
    """

    if [ "${optimizeBins_method}" = 'PBO' ] || [ "${optimizeBins_method}" = 'DBO' ]; then

        cp -rf ${bin_floder} ./${id}_Best_Bin
        cp -rf *${bin_qs} ./${id}_Best_Bin_quality_report.tsv

    else

        pick_best_bin4PBO_DBO.py -pb ${id}_PBO_best_bins -ps ${id}_PBO_quality_report.tsv -bd ${id}_DASTool_bins  -bs quality_report.tsv -s ${id} -a ${params.completeness} -b ${params.contamination} -r ${params.similarity_ratio}

    fi

    """

}