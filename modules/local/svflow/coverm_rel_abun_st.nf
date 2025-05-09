
process COVERMRELABUNST {

    label 'process_medium'

    input:
    path(bins_rel_abun)

    output:
    path("${params.pipeline_prefix}_CoverM_bins_st_rel_abun.xls"),emit:"st_rel_abun"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    SV_coverm_result_normalization.py ${bins_rel_abun} ${params.pipeline_prefix}_CoverM_bins_st_rel_abun.xls
    

    """

}
