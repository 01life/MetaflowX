
process ANTISMASH {

    tag "$id"

    label 'process_single'

    input:
    tuple val(id),path(contigs)

    output:
    tuple val(id),path("all.bin",type:'dir'),emit:"all_bin"
    path("${id}.zip"),emit:"res"

    when:
    task.ext.when == null || task.ext.when

    script:
    def antismash_options = params.antismash_options ?: ""
    """

    ln -s ${contigs} ./all.bin.fa
    antismash all.bin.fa -c ${task.cpus} ${antismash_options}

    mv all.bin/all.bin.zip ${id}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash-lite: \$(echo \$(antismash --version) | sed 's/antiSMASH //')
    END_VERSIONS    
    
    """

    stub:
    """
    mkdir all.bin
    touch ${id}.zip

    """

}
