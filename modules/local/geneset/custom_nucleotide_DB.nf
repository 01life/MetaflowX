
process CUSTOMNTDB {

    label 'process_single'
    
    input:
    path(cds)
    path(abundance)
    path(nucleotide_db)

    output:
    path("${params.pipeline_prefix}_geneset_function_${params.ntDB_name}*.xls")
    path("customnt_report"),emit:"customnt_info"

    when:
    task.ext.when == null || task.ext.when

    script: 
    def bowtie2_options = params.ntDB_bowtie2_options ?:  ""
    
    """

    bowtie2-build --threads ${task.cpus} ${nucleotide_db} nt.index

    bowtie2 -f -x nt.index -U ${cds} -p ${task.cpus} ${bowtie2_options} | get_align_source4Sam.py > ${params.pipeline_prefix}_geneset_function_${params.ntDB_name}_annotation.xls

    if [ -s "${params.pipeline_prefix}_geneset_function_${params.ntDB_name}_annotation.xls" ]; then

        functional_profile.pl -k ${params.pipeline_prefix}_geneset_function_${params.ntDB_name}_annotation.xls -f ${abundance} -o ./ -pr ${params.pipeline_prefix}_geneset_function_${params.ntDB_name}_abundance.xls

        mkdir customnt_report
        cp ${params.pipeline_prefix}_geneset_function_${params.ntDB_name}_abundance.xls customnt_report/customNT.xls

    else
        touch ${params.pipeline_prefix}_geneset_function_${params.ntDB_name}_NO_RESULT_abundance.xls
        touch ${params.pipeline_prefix}_geneset_function_${params.ntDB_name}_NO_RESULT_annotation.xls
        touch customnt.report.list
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS

    """
}
