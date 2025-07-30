
process CONTIGFILTER {
    
    tag "$id"

    label 'process_low'
    
    input:
    tuple val(id),path(contig,stageAs:'tmp/*')

    output:
    tuple val(id),path("${id}_contig_${params.min_contig_len}*.fa"),emit:"contigs"
    tuple val(id),path("${id}_contig_rename_mapping.xls"),emit:"contigs_map"

    when:
    task.ext.when == null || task.ext.when

    script:
   
    """

    fasta=\$(ls tmp/* 2>/dev/null | head -n 1)
    if [ -z "\$fasta" ]; then
        echo "Error: No ${id} contig files found in tmp/ directory"
        exit 1
    fi

    seqtk seq -L ${params.min_contig_len} \$fasta > seq.contigs.fa

    contig_rename.py -t seq.contigs.fa -s ${id} -m ${id}_contig_rename_mapping.xls -o ${id}_contig_${params.min_contig_len}.fa

    """
    stub:
    """
    touch ${id}_contig_${params.min_contig_len}_1.fa
    touch ${id}_contig_rename_mapping.xls
    """ 
}
