process CATEVALBYID {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(eval)

    output:
    tuple val(id),path("${id}_allBins.eval.txt"), emit:eval

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${eval} > ${id}_allBins.eval.txt

    """
    stub:
    """
    touch ${id}_allBins.eval.txt
    """

}