process  PICKPBODBO {
    
    tag "$id"
    
    label 'process_low'

    input:
    val(optimizeBins_method)
    tuple val(id), path(bin_floder),path(bin_qs)

    output:
    tuple val(id),path("${id}/${id}_Best_Bin"),path("${id}/${id}_Best_Bin_quality_report.tsv"),emit:"best", optional: true
    tuple val(id),path("${id}/${id}_Best_Bin"),emit:"bestBin", optional: true
    tuple val(id),path("${id}/${id}_pick_PBO_DBO_error.log"),emit:"errorlog", optional: true
    
    when:
    task.ext.when == null || task.ext.when

    script:
    
    """
    mkdir ${id}
    
    if [ "${optimizeBins_method}" = 'PBO' ] || [ "${optimizeBins_method}" = 'DBO' ]; then
        # check file and floder exit
        if [ -d "${bin_floder}" ] ; then
            cp -rf ${bin_floder} ${id}/${id}_Best_Bin
            cp -rf *${bin_qs} ${id}/${id}_Best_Bin_quality_report.tsv
        else
            echo "Not PBO or DBO Required directories or files do not exist." > ${id}/${id}_pick_PBO_DBO_error.log
        fi
    else
        # check file and floder exit
        if [ -d "${id}_PBO_best_bins" ] && [ -f "${id}_PBO_quality_report.tsv" ] && [ -d "${id}_DASTool_bins" ] && [ -f "quality_report.tsv" ]; then
            pick_best_bin4PBO_DBO.py -pb ${id}_PBO_best_bins -ps ${id}_PBO_quality_report.tsv -bd ${id}_DASTool_bins -bs quality_report.tsv -s ${id} -a ${params.completeness} -b ${params.contamination} -r ${params.similarity_ratio}
            mv ${id}_Best_Bin ${id}_Best_Bin_quality_report.tsv ${id}
        else
            echo "Not PBO and DBO  Required directories or files do not exist." > ${id}/${id}_pick_PBO_DBO_error.log
        fi
    fi

    """

}