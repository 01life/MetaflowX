
process CUSTOMPRODB {

    label 'process_single'

    input:
    path(pep)
    path(abundance)
    path(protein_db)

    output:
    path("${params.pipeline_prefix}_geneset_function_${params.proDB_name}*.xls")
    path("custompr_report"),emit:"custompr_info"

    when:
    task.ext.when == null || task.ext.when

    script:
    def proDB_diamond_options = params.proDB_diamond_options ?: "" 
    
    """

    diamond makedb --in ${protein_db} -d index

    diamond blastp -d index -k 1 -q ${pep} -o ${params.proDB_name}.diamond.out ${proDB_diamond_options} -p ${task.cpus}

    cut -f 1,2 ${params.proDB_name}.diamond.out > ${params.pipeline_prefix}_geneset_function_${params.proDB_name}_annotation.xls

    if [ -s "${params.pipeline_prefix}_geneset_function_${params.proDB_name}_annotation.xls" ]; then

	    functional_profile.pl -k ${params.pipeline_prefix}_geneset_function_${params.proDB_name}_annotation.xls -f ${abundance} -o ./ -pr ${params.pipeline_prefix}_geneset_function_${params.proDB_name}_abundance.xls
        
        mkdir custompr_report
        cp ${params.pipeline_prefix}_geneset_function_${params.proDB_name}_abundance.xls custompr_report/customPR.xls
     
    else
        touch ${params.pipeline_prefix}_geneset_function_${params.proDB_name}_NO_RESULT_abundance.xls
        touch ${params.pipeline_prefix}_geneset_function_${params.proDB_name}_NO_RESULT_abundance.xls
        touch custompr.report.list 
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(echo \$(diamond --version 2>&1) | sed -e 's/diamond version //g')
    END_VERSIONS

    """
}
