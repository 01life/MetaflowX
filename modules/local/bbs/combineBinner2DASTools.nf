process COMBINEBINNER {
    
    label 'process_low'

    input:
    path(contig2bin)
    path(protein_list)
    path(contig_list)

    output:
    path("binner_combination.txt"),emit:"binner_combination"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    combineBinner2DASTools.py -t ${contig2bin} -p ${protein_list} -c ${contig_list}

    """
    
}