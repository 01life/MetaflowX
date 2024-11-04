
process CONTIGFILTER {
    
    tag "$id"

    label 'process_low'
    
    input:
    tuple val(id),path(contig)

    output:
    tuple val(id),path("${id}_contig_${params.min_contig_len}*.fa"),emit:"contigs"
    path("${id}_contig_rename_mapping.xls")

    when:
    task.ext.when == null || task.ext.when

    script:
   
    """
        
    seqtk seq -L ${params.min_contig_len} ${contig} > seq.contigs.fa
    
    if [ -f "${id}_contig_${params.min_contig_len}.fa" ]; then
        timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')  
        contig_rename.py -t seq.contigs.fa -s ${id} -m ${id}_contig_rename_mapping.xls -o ${id}_contig_${params.min_contig_len}_\$timestamp.fa
    else
        contig_rename.py -t seq.contigs.fa -s ${id} -m ${id}_contig_rename_mapping.xls -o ${id}_contig_${params.min_contig_len}.fa
    fi

    """
}
