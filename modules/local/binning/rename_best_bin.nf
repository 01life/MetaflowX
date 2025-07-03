process  RENAMEBEXTBIN {
    
    tag "$id"
    
    label 'process_low'

    input:
    val(optimizeBins_method)
    tuple val(id), path(bin_floder),path(bin_qs)

    output:
    tuple val(id),path("${id}_${optimizeBins_method}_Best_Bin"),emit:"bestBin", optional: true
    tuple val(id),path("${id}_${optimizeBins_method}_Best_Bin_quality_report.tsv"),emit:"bestBinQS", optional: true
    
    when:
    task.ext.when == null || task.ext.when

    script:
    
    """
    cp -rf ${bin_floder} ${id}_${optimizeBins_method}_Best_Bin
    cp -rf ${bin_qs} ${id}_${optimizeBins_method}_Best_Bin_quality_report.tsv

    """

}