
process MERGEBGCPROFILE {
    
    label 'process_single'
    
    input:
    path(bgcfamily)
    path(samplesheet)

    output:
    path("${params.pipeline_prefix}_*.xls"),optional:true
    path("bigmap_report"),emit:"bigmap_info",optional:true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ls *.xls |awk -F "_BiG-MAP" '{print\$1"\t./"\$0}' > all.BiG-MAP

    parallel "grep {1} all.BiG-MAP > {1}.list ; merge_multi_matrix_XY_0.py {1}.list ${params.pipeline_prefix}_{1}_raw.xls; col_reorder.pl ${params.pipeline_prefix}_{1}_raw.xls 1 1 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_{1}.xls" ::: BiG-MAP_corecov BiG-MAP_coreRAW BiG-MAP_coreRPKM BiG-MAP_coreTPM BiG-MAP_RAW BiG-MAP_RPKM BiG-MAP_TPM BiG-MAP_cov

    # Delete result files before sorting.
    rm -f ${params.pipeline_prefix}_*_raw.xls

    mkdir bigmap_report
    cp ${params.pipeline_prefix}_BiG-MAP_RPKM.xls bigmap_report/antismash.xls
    
    """
    
}
