
process CDHIT {

    label 'process_single'

    input:
    val(sample_number)
    val(finish_number)
    path(cds)
    path(faa)

    output:
    path("${params.pipeline_prefix}_geneset_protein.fa"),emit:"pep"
    path("${params.pipeline_prefix}_geneset_gene.fa"),emit:"cds"
    path("${params.pipeline_prefix}_geneset_gene_length.xls"),emit:"gene_length"
    path("${params.pipeline_prefix}_geneset_gene_info.xls"),emit:"gene_info"
    path("${params.pipeline_prefix}_geneset_cdhit_clstr.txt"),emit:"clstr"
    path("split/*.fa"),emit:"split"
    path("geneset_Gene_report"),emit:"geneset_gene_report"

    when:
    sample_number == finish_number

    script:
    def options = params.cdhit_options ?: ""
    def cdhitOptions = "\" ${options} \""
    def queueOptions = "\" -A ${params.Account} -p ${task.queue} -c ${params.cdhit_split_run_eachthread} \""
    
    """

    cat ${cds} > all_prodigal.cds
    cat ${faa} > all_prodigal.faa
    
    gene_stat.py all_prodigal.cds geneset_gene_stat.txt

    mkdir geneset_Gene_report
    cp geneset_gene_stat.txt geneset_Gene_report/genesetSampleStat.txt

    filter_GeneLength.py -i all_prodigal.cds -o all.L-150.fa -l ${params.gene_min_length} -a ${cdhitOptions} -b ${queueOptions} -t ${task.cpus} -x ${params.cdhit_split_run_eachthread} -s ${params.cdhit_split_run_threshold} -c ${params.cdhit_geneset_chunk_size} -T ${task.executor}

    sh cd-hit.sh
    

    mv unique.fa.clstr ${params.pipeline_prefix}_geneset_cdhit_clstr.txt

    rename_GeneSet.py -n unique.fa -p all_prodigal.faa -q ${params.pipeline_prefix}

    seqkit split ${params.pipeline_prefix}_geneset_protein.fa -s ${params.eggnog_protein_chunk_size} -O split

    cp ${params.pipeline_prefix}_geneset_gene_length.xls geneset_Gene_report/genesetLenStat.xls
    
    rm -rf all*


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit-est 2>&1 | grep 'CD-HIT version' | sed 's/CD-HIT //')
    END_VERSIONS

    """
    stub:
    """
    mkdir split
    touch split/gene1.fa
    touch ${params.pipeline_prefix}_geneset_protein.fa
    touch ${params.pipeline_prefix}_geneset_gene.fa
    touch ${params.pipeline_prefix}_geneset_gene_info.xls
    touch ${params.pipeline_prefix}_geneset_cdhit_clstr.txt
    touch ${params.pipeline_prefix}_geneset_gene_length.xls
    mkdir geneset_Gene_report
    """

}
