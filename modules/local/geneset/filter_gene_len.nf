
process GENEFILTER {

    label 'process_single'

    input:
    val(sample_number)
    val(finish_number)
    path(cds)
    path(faa)

    output:
    path("*.multi.cdhit.task"),emit:"multi_task" , optional: true
    path("single.cdhit.task.num"),emit:"single_task" , optional: true
    path("all.L-*.fa"),emit:"allcds"
    path("all_prodigal.faa"),emit:"allpep"
    path("geneset_Gene_report"),emit:"geneset_gene_report"


    when:
    sample_number == finish_number

    script:
    
    """
    cat ${cds} > all_prodigal.cds
    cat ${faa} > all_prodigal.faa
    
    gene_stat.py all_prodigal.cds geneset_gene_stat.txt

    mkdir geneset_Gene_report
    cp geneset_gene_stat.txt geneset_Gene_report/genesetSampleStat.txt

    filter_GeneLength.py -i all_prodigal.cds -o all.L-150.fa -l ${params.gene_min_length}  -s ${params.cdhit_split_run_threshold} -c ${params.cdhit_geneset_chunk_size} 


    """

}
