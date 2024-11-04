process COMBINEBINNER {
    
    tag "$id"
    
    label 'process_single'

    input:
    tuple val(id),path(contigs2tsv),path(depth),path(contig)

    output:
    tuple val(id),path("*__*tsv"),emit:"combine_bin_tsv"
    tuple val(id),path("*:*",type:'dir'),emit:"combine_bin_fa"
    tuple val(id),path("*allcontigs2bin.txt"),emit:"allcontigs2bin"

    when:
    task.ext.when == null || task.ext.when

    script:
    
    """
    comebine_bins_refine.py -info ${depth} -bins ${contigs2tsv} -t ${task.cpus} -i ${id} -c ${contig}
    
    """

}