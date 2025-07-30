
process SINGLECDHIT {

    label 'process_single'

    input:
    path(allcds)
    path(allpep)
    path(task_num)

    output:
    path("*_geneset_protein.fa"),emit:"pep"
    path("*_geneset_gene.fa"),emit:"cds"
    path("*_geneset_gene_length.xls"),emit:"gene_length"
    path("*_geneset_gene_ID_mapping.xls"),emit:"gene_info"
    path("*_geneset_cdhit_clstr.txt"),emit:"clstr"
    path("split/*.fa"),emit:"split"
    path("geneset_Gene_report"),emit:"geneset_gene_report"

    script:
    def options = params.cdhit_options ?: ""
    
    """

    cat ${allcds} > all.L-150.fa

    cd-hit-est -T ${task.cpus} -i all.L-150.fa -o unique.fa ${options}

    rm -rf all.L-150.fa

    mv unique.fa.clstr ${params.pipeline_prefix}_geneset_cdhit_clstr.txt

    merge_psort_pep.pl -fa unique.fa -input ${allpep} \\
        -preprocessed -cpu ${task.cpus} \\
        -output unique_protein.fa

    awk '/^>/ { printf(">NC_%010d\\n", ++count); next } { print }' unique.fa > ${params.pipeline_prefix}_geneset_gene.fa
    paste <(grep ">" unique.fa | cut -d '>' -f2 | cut -d ' ' -f1) <(grep ">" ${params.pipeline_prefix}_geneset_gene.fa | cut -d '>' -f2) > ${params.pipeline_prefix}_geneset_gene_ID_mapping.xls

    awk '/^>/ { printf(">NC_%010d\\n", ++count); next } { print }' unique_protein.fa > ${params.pipeline_prefix}_geneset_protein.fa

    rm -rf unique.fa unique_protein.fa

    seqkit fx2tab -l -i -n ${params.pipeline_prefix}_geneset_gene.fa | awk -v OFS='\\t' '{print NR, \$0}' | sed '1iID\tname\tlength' > ${params.pipeline_prefix}_geneset_gene_length.xls
    
    seqkit split ${params.pipeline_prefix}_geneset_protein.fa -s ${params.eggnog_protein_chunk_size} -O split

    mkdir geneset_Gene_report
    cp ${params.pipeline_prefix}_geneset_gene_length.xls geneset_Gene_report/genesetLenStat.xls
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit-est 2>&1 | grep 'CD-HIT version' | sed 's/CD-HIT //')
    END_VERSIONS

    """

    stub:
    """
    touch ${params.pipeline_prefix}_geneset_protein.fa
    touch ${params.pipeline_prefix}_geneset_gene.fa
    touch ${params.pipeline_prefix}_geneset_gene_length.xls
    touch ${params.pipeline_prefix}_geneset_cdhit_clstr.txt
    touch ${params.pipeline_prefix}_geneset_gene_ID_mapping.xls
    mkdir split
    mkdir geneset_Gene_report
    touch split/chunk1.fa
    """

}
