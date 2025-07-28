
process GENEPROFILE {

    label 'process_low'

    input:
    path(annotation)
    path(merge)
    path(cog_db_category)          
    path(go_db_category)
    path(kegg_db_category)           
    path(cazy_db_category)

    output:
    path("${params.pipeline_prefix}_geneset_function_emapper_org_annotation.xls"),emit:"All_anontations"
    path("${params.pipeline_prefix}_geneset_gene_abundance.xls"),emit:"total_abundance"
    path("genesetFunction_report"),emit:"geneset_function_report"
    path("${params.pipeline_prefix}_geneset*.xls")
    path("${params.pipeline_prefix}_geneset*abundance.xls"),emit:"func_abun"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    cat ${annotation} > ${params.pipeline_prefix}_geneset_function_emapper_org_annotation.xls
    
    get_abundance_4eggNOG-mapper.py -e ${params.pipeline_prefix}_geneset_function_emapper_org_annotation.xls -f ${merge} -q ${params.pipeline_prefix}


    mkdir genesetFunction_report

    # check COG file exist or not
    if [ -f "${params.pipeline_prefix}_geneset_function_cog_catF_abundance.xls" ]; then
        cogLevel1_stat.py ${params.pipeline_prefix}_geneset_function_cog_catF_abundance.xls ${cog_db_category}  cog_level1_sort.txt
        cp cog_level1_sort.txt genesetFunction_report/genesetCOG.txt
    fi

    # check GO file exist or not
    if [ -f "${params.pipeline_prefix}_geneset_function_GOs_abundance.xls" ]; then
        goLevel1_stat.py ${params.pipeline_prefix}_geneset_function_GOs_abundance.xls ${go_db_category} go_level1_sort.txt
        cp go_level1_sort.txt genesetFunction_report/genesetGO.txt
    fi

    # check KEGG file exist or not
    if [ -f "${params.pipeline_prefix}_geneset_function_KEGG_ko_abundance.xls" ]; then
        keggLevel1_stat.py ${params.pipeline_prefix}_geneset_function_KEGG_ko_abundance.xls ${kegg_db_category} kegg_level1_sort.txt
        cp kegg_level1_sort.txt genesetFunction_report/genesetKEGG.txt
    fi

    # check CAZY file exist or not
    if [ -f "${params.pipeline_prefix}_geneset_function_CAZy_abundance.xls" ]; then
        cazyLevel1_stat.py ${params.pipeline_prefix}_geneset_function_CAZy_abundance.xls ${cazy_db_category} cazy_level1_sort.txt
        cp cazy_level1_sort.txt genesetFunction_report/genesetCAZY.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS

    """


    stub:
    """
    touch ${params.pipeline_prefix}_geneset_function_emapper_org_annotation.xls
    touch ${params.pipeline_prefix}_geneset_gene_abundance.xls
    mkdir genesetFunction_report
    touch genesetFunction_report/genesetCOG.txt
    touch genesetFunction_report/genesetGO.txt
    touch genesetFunction_report/genesetKEGG.txt
    touch genesetFunction_report/genesetCAZY.txt
    touch ${params.pipeline_prefix}_geneset_COG_abundance.xls
    touch ${params.pipeline_prefix}_geneset_GO_abundance.xls
    touch ${params.pipeline_prefix}_geneset_KEGG_abundance.xls
    touch ${params.pipeline_prefix}_geneset_CAZy_abundance.xls
    touch ${params.pipeline_prefix}_geneset_function_cog_catF_abundance.xls
    touch ${params.pipeline_prefix}_geneset_function_GOs_abundance.xls
    touch ${params.pipeline_prefix}_geneset_function_KEGG_ko_abundance.xls
    touch ${params.pipeline_prefix}_geneset_function_CAZy_abundance.xls
    touch ${params.pipeline_prefix}_geneset_function_cog_category.xls
    touch ${params.pipeline_prefix}_geneset_function_GOs_category.xls
    touch ${params.pipeline_prefix}_geneset_function_KEGG_ko_category.xls       
    touch ${params.pipeline_prefix}_geneset_function_CAZy_category.xls
    """
}
