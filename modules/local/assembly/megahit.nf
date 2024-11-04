process MEGAHIT {
    
    tag "$id"

    label 'process_single'

    input:
    tuple val(id),path(reads)

    output:
    tuple val(id),path("${id}/${id}_contigs.fa"),emit:"contigs"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.megahit_options ?: ""
    def illumina_reads = reads ? ( params.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}" ) : ""
    
    """
    megahit ${illumina_reads} -o ${id} -t ${task.cpus} ${options}
    seqtk rename ${id}/final.contigs.fa "${id}|" > ${id}/${id}_contigs.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS    

    """
    
    stub:
    """
    mkdir ${id}
    touch ${id}/${id}_contigs.fa

    """

}
