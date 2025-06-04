
process METASPADESM128 {
    
    tag "$id"

    input:
    tuple val(id),path(reads)

    output:
    tuple val(id),path("${id}_contigs.fa"), emit: "contigs", optional: true
    tuple val(id),path("${id}_scaffolds.fa"), emit: "scaffolds", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.metaspades_options ?: ""
    def maxmem = task.memory.toGiga()
    def illumina_reads = reads ? ( params.single_end ? "-s ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}" ) : ""
    
    """
    
    spades.py ${illumina_reads} -o ${id} -t ${task.cpus} --memory $maxmem ${options} || rm -rf ${id} 
    
    if [ -e "${id}/contigs.fasta" ]; then
    
        seqtk seq -L ${params.min_contig_len} ${id}/contigs.fasta > ${id}/seq.contigs.fa
    
        seqtk rename ${id}/seq.contigs.fa "${id}|" > ${id}_contigs.fa

        if [ -e "${id}/scaffolds.fasta" ]; then

            seqtk seq -L ${params.min_contig_len} ${id}/scaffolds.fasta > ${id}/seq.scaffolds.fasta
            seqtk rename ${id}/seq.scaffolds.fasta "${id}|" > ${id}_scaffolds.fa
        
        fi

        rm -rf ${id} 
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaspades: \$(metaspades.py --version | sed "s/SPAdes genome assembler v//; s/ \\[.*//")
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS

    """

    stub:
    """
    touch ${id}_contigs.fa

    """
}
