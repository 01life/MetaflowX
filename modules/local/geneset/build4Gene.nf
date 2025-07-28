
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

    stub:
    """
    mkdir index
    touch index/index.1.bt2
    touch index/index.2.bt2
    touch index/index.3.bt2
    touch index/index.4.bt2
    touch index/index.rev.1.bt2
    touch index/index.rev.2.bt2
    """

}
