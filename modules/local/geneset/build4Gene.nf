
process BUILD4GENE {

    label 'process_single'

    input:
    path(cds)

    output:
    path("index*"),emit:"index"

    when:
    task.ext.when == null || task.ext.when

    script:
    
    """    
    bowtie2-build --threads ${task.cpus} ${cds} index

    """

}
