
process VFDB {

    label 'process_single'

    input:
    path(pep)
    path(abundance)
    path(VFDB_db)

    output:
    path("${params.pipeline_prefix}_geneset_function_VFDB*.xls")
    path("${params.pipeline_prefix}_geneset_function_VFDB_annotation.xls"),emit:"vfdb_anno"
    path("vfdb_report"),emit:"vfdb_info"

    when:
    task.ext.when == null || task.ext.when

    script:
    def proDB_diamond_options = params.proDB_diamond_options ?: "" 
    
    """

    diamond makedb --in ${VFDB_db} -d index

    diamond blastp -d index -k 1 -q ${pep} -o VFDB.diamond.out ${proDB_diamond_options} -p ${task.cpus}

    cut -f 1,2 VFDB.diamond.out > ${params.pipeline_prefix}_geneset_function_VFDB_annotation.xls

	functional_profile.pl -k ${params.pipeline_prefix}_geneset_function_VFDB_annotation.xls -f ${abundance} -o ./ -pr ${params.pipeline_prefix}_geneset_function_VFDB_abundance.xls

    mkdir vfdb_report
    cp ${params.pipeline_prefix}_geneset_function_VFDB_abundance.xls vfdb_report/VFDB.xls

    #echo VFDB \$PWD/${params.pipeline_prefix}_geneset_function_VFDB_abundance.xls > vfdb.report.list

    """

    stub:
    """
    touch ${params.pipeline_prefix}_geneset_function_VFDB_annotation.xls
    mkdir vfdb_report
    touch vfdb_report/VFDB.xls
    touch vfdb.report.list
    """
}
