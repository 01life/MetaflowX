
process CDHIT {

    label 'process_low'

    input:

    path(orig_pep)
    path(orig_cds)
    path(orig_gene_length)
    path(orig_gene_info)
    path(orig_clstr)
    path(orig_split)
    path(orig_report)


    output:
    path("${params.pipeline_prefix}_geneset_protein.fa"),emit:"pep"
    path("${params.pipeline_prefix}_geneset_gene.fa"),emit:"cds"
    path("${params.pipeline_prefix}_geneset_gene_length.xls"),emit:"gene_length"
    path("${params.pipeline_prefix}_geneset_gene_info.xls"),emit:"gene_info"
    path("${params.pipeline_prefix}_geneset_cdhit_clstr.txt"),emit:"clstr"
    path("split/*.fa"),emit:"split"
    path("geneset_Gene_report"),emit:"geneset_gene_report"



    script:
    
    """
    cp ${orig_pep} ${params.pipeline_prefix}_geneset_protein.fa

    cp ${orig_cds} ${params.pipeline_prefix}_geneset_gene.fa
    
    cp ${orig_gene_length} ${params.pipeline_prefix}_geneset_gene_length.xls
    
    cp ${orig_gene_info} ${params.pipeline_prefix}_geneset_gene_info.xls
    
    cp ${orig_clstr} ${params.pipeline_prefix}_geneset_cdhit_clstr.txt

    mkdir split
    cp ${orig_split} ./split
    
    cp -rf ${orig_report} geneset_Gene_report

    """

}
