
process BOWTIE2GENE {
    
    tag "$id"

    label 'process_single'

    input:
    path(index)
    path(gene_length)
    tuple val(id),path(reads)

    output:
    path("${id}_abundance.xls"),emit:"abundance"
    path("${id}_total_abundance.xls"),emit:"total_abun"
    path("${id}_geneset_bowtie2_log.txt"),emit:"bowtie2_log"

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_args = params.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def options = params.geneset_profile_bowtie2_options ?: ""

    """

    bowtie2 -x index ${reads_args} -S ${id}.sam -p ${task.cpus} ${options} > ${id}_geneset_bowtie2_log.txt 2>&1

    # Notice: please not use parameter -nz to run get_geneset_abundance.py!
    get_geneset_abundance.py -i ${id}.sam -g ${gene_length} -n ${task.cpus} -ng -p ${id} -s ${id}

    #Clean up intermediate files.
    rm -rf *.sam
    #rm \$(readlink -f *.bt2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS

    """

    stub:
    """
    touch ${id}_abundance.xls
    touch ${id}_geneset_bowtie2_log.txt
    
    """

}
