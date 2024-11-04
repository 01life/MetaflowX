
process PASTERESULT {


    label 'process_single'

    input:

    path(subclstr)
    path(uniq_fa)
    path(allpep)

    output:
    path("*_geneset_protein.fa"),emit:"pep"
    path("*_geneset_gene.fa"),emit:"cds"
    path("*_geneset_gene_length.xls"),emit:"gene_length"
    path("*_geneset_gene_info.xls"),emit:"gene_info"
    path("*_geneset_cdhit_clstr.txt"),emit:"clstr"
    path("*split/*.fa"),emit:"split"
    path("*geneset_Gene_report"),emit:"geneset_gene_report"

    script:


    """
    prefix='multi_cdhit'

    cat ${subclstr} > unique.fa.clstr

    mv unique.fa.clstr \$prefix\\_geneset_cdhit_clstr.txt

    rename_GeneSet.py -n  ${uniq_fa} -p ${allpep} -q \$prefix

    seqkit split \$prefix\\_geneset_protein.fa -s ${params.eggnog_protein_chunk_size} -O \$prefix\\_split

    mkdir \$prefix\\_geneset_Gene_report
    cp \$prefix\\_geneset_gene_length.xls \$prefix\\_geneset_Gene_report/genesetLenStat.xls
    
    """

}
