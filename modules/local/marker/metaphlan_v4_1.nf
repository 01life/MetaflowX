
process METAPHLANV41 {
    
    tag "$id"
    
    label 'process_high'

    input:
    tuple val(id),path(reads)
    path(mpa_db)

    output:
    tuple val(id),path("${id}.xls"),emit:"profile"
    tuple val(id),path("${id}_mpa_bowtie2.bz2"),emit:"bowtie2out"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.mpa_options ?: ""
    def input = params.single_end ? "${reads}" : "${reads[0]},${reads[1]}"
    """
    metaphlan ${input} --bowtie2out ${id}_mpa_bowtie2.bz2 --output_file ${id}.xls --nproc ${task.cpus} --bowtie2db ${mpa_db} --index ${params.mpa_index} --sample_id ${id} ${options}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
        
    """
    
    stub:
    """
    echo "Metaphlanv41 not run, stub test" > ${id}.xls 
    touch ${id}_mpa_bowtie2.bz2

    """
}
