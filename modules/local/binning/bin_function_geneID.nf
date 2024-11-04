
process BINFUNCTIONGENEID {

    label 'process_low'

    input:
    path(gene_info)
    path(clstr)
    path(annotations)
    path(genomes)
    path(cog_db_category)          
    path(go_db_category)
    path(kegg_db_category)           
    path(cazy_db_category)

    output:
    path("${params.pipeline_prefix}*.xls"),emit:"functional" , optional: true
    path("eachBinFunction_report"),emit:"eachBinFunction" , optional: true
    path("get_bin_function_error_log.txt"),emit:"warning" , optional: true
    

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    get_eachBin_FunctionGeneID.py -e ${annotations} \
    -i ${gene_info} \
    -c ${clstr}  \
    -b  ./ \
    -q  ${params.pipeline_prefix}_bin_function \
    -o  ./

    paste_eachBinFunction2animation.py -k ${params.pipeline_prefix}_bin_function_KEGG_ko_annotation.xls -K ${kegg_db_category} -g ${params.pipeline_prefix}_bin_function_GOs_annotation.xls -G ${go_db_category}  -c ${params.pipeline_prefix}_bin_function_cog_catF_annotation.xls -C ${cog_db_category} -a ${params.pipeline_prefix}_bin_function_CAZy_annotation.xls -A ${cazy_db_category}
    

    if [ \$(wc -l < eachBin_Function_stat.xls) -lt 2 ]; then
        
        timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')
    
        cat <<-OUTLOG > get_bin_function_error_log.txt
    
==========Start at : `date` ==========
### Step ${task.process}
Oppo~ Do not get any annotation of all bins. Please check your input and data. SO BAD~~~ T -_- T
==========End at : `date` ==========

OUTLOG

    else

        mkdir eachBinFunction_report
        cp eachBin_Function_stat.xls eachBinFunction_report/eachBinFunction.xls

    fi


    #echo eachBinFunction \$PWD/eachBin_Function_stat.xls > eachBinFunction.report.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS

    """


}
