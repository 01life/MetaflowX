
process BINSPECIFIEDFUNCTION {

    label 'process_low'

    input:
    path(gene_info)
    path(clstr)
    path(annotations)
    path(genomes)
    val(database)

    output:
    path("${params.pipeline_prefix}_bin_function_${database}_annotation.xls"),emit:"functional" , optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    bf_get_eachBin_specified_FunctionGeneID.py -e ${annotations} \
    -i ${gene_info} \
    -c ${clstr}  \
    -b  ./ \
    -q  ${params.pipeline_prefix}_bin_function \
    -d  ${database} \
    -o  ./



    """


}
