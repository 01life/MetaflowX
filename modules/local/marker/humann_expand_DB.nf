
process HUMANNEXPAND {
    
    tag "$id"
    
    label 'process_medium'

    input:
    tuple val(name),val(expand_db),val(id),path(profile)

    output:
    tuple val(name),path("${id}_${name}_relab.xls"),emit:"profile"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.humann_regroup_table_options ?: ""
    """

    humann_regroup_table ${options} -i ${profile} -c ${expand_db} -o ${id}_${name}_relab.xls
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$( humann --version | sed 's/humann //' )
    END_VERSIONS
        
    """

}
