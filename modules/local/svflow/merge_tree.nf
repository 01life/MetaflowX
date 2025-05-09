
process MERGETREE {
    
    label 'process_single'
    
    input:
    path(tree)
    path(bins_taxonomic)
    path(bins_info)

    output:
    path("${params.pipeline_prefix}_tree_info.xls"),emit:"tree"

    when:
    task.ext.when == null || task.ext.when

    script:
    
    """
    
    SV_get_mergetree.py -i ${tree} -t ${bins_taxonomic} -a ${bins_info} -o ${params.pipeline_prefix}_tree_info.xls

    
    """
    
}
