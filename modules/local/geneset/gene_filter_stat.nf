
process GENEFILTERSTAT {

    tag "$id"

    label 'process_single'

    input:
    tuple val(id),path(cds),path(faa)

    output:
    path("${id}_L${params.gene_min_length}_cds.fa"),emit:"cds"
    path("${id}_L${params.gene_min_length}_protein.fa"),emit:"pep"
    path("${id}_L${params.gene_min_length}_gene_stat.xls"),emit:"gene_stat"
    path("${id}_L${params.gene_min_length}_gene_info.xls"),emit:"gene_info"
    path("${id}_L${params.gene_min_length}_gene_length.xls"),emit:"gene_length"

    script:
    """
    gene_stat.py -q ${id} -n ${cds} -p ${faa} -l ${params.gene_min_length}

    """

}
