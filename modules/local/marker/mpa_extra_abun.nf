
process MPAEXTRAABUN {
    
    tag "$id-$method"
    
    label 'process_single'

    input:
    tuple val(method),val(id),path(bowtie2out)
    path(mpa_db)

    output:
    tuple val(method),path("${id}_${method}.xls"),emit:"profile"
    tuple val(id),path("${id}_${method}.xls"),emit:"out4humann"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    metaphlan ${bowtie2out} -t ${method} --input_type bowtie2out --bowtie2db ${mpa_db} --index ${params.mpa_index} --nproc ${task.cpus} -o ${id}_${method}.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
        
    """

}
