
process PRODIGAL {
    
    tag "$id"
    
    label 'process_single','error_ignore'

    input:
    tuple val(id),path(contigs)

    output:
    tuple val(id),path("${id}_gene.fa"),emit:"cds"
    tuple val(id),path("${id}_protein.fa"),emit:"faa"
    path("${id}_gene.fa"),emit:"noid_cds"
    path("${id}_protein.fa"),emit:"noid_faa"


    when:
    contigs.size() > 0

    script:
    def options = params.prodigal_options ?: ""
    
    """
    prodigal -i ${contigs} -o ${id}_gene.coords.gbk -a ${id}_protein.fa -d ${id}_gene.fa ${options} 2>prodigal.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
    END_VERSIONS

    """

    stub:
    """
    echo ">${id}_gene" > ${id}_gene.fa
    echo "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG" >> ${id}_gene.fa
    echo ">${id}_protein" > ${id}_protein.fa
    echo "MDGSFCTGAGS" >> ${id}_protein.fa
    """
}
