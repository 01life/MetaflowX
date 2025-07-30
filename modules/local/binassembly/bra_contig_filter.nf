
process BRACONTIGFILTER {
    
    tag "$id"

    label 'process_low'
    
    input:
    tuple val(id),path(contig)

    output:
    path("${id}_reassembly_contigs_${params.min_contig_len}.fa"),emit:"filtercontigs"
    path("${id}_reassembly_contigs_rename_mapping.xls"),emit:"renamecontigs"

    when:
    task.ext.when == null || task.ext.when

    script:
   
    """
        
    bra_spades_length_cov_filter.py -f ${contig}  -l ${params.reassembly_min_contig_len}   -p reassembly_bin

    #_remove_outlier.fa

    contig_rename.py -t reassembly_bin_remove_outlier.fa -s ${id} -m ${id}_reassembly_contigs_rename_mapping.xls -o  ${id}_reassembly_contigs_${params.min_contig_len}.fa

    """
    stub:
    """
    touch ${id}_reassembly_contigs_${params.min_contig_len}.fa
    touch ${id}_reassembly_contigs_rename_mapping.xls
    """

}
