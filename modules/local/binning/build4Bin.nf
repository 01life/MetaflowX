
process BUILD4BIN {

    label 'process_medium'

    input:
    path(genomes)

    output:
    path("all.bin.fa*"),emit:"bowtie2Index"

    when:
    task.ext.when == null || task.ext.when

    script:
    
    """    

    cat ${genomes} > all.bin.fa

    # bowtie2
    bowtie2-build -f all.bin.fa all.bin.fa --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS

    """
    
    stub:
    """
    touch all.bin.fa
    touch all.bin.fa.1.bt2
    touch all.bin.fa.2.bt2
    touch all.bin.fa.3.bt2
    touch all.bin.fa.4.bt2
    touch all.bin.fa.rev.1.bt2
    touch all.bin.fa.rev.2.bt2
    """

}
