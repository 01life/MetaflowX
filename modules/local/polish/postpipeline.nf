process POSTPIPELINE {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(reads),path(report)

    output:
    path("${id}/*")
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """     
    
    mkdir ${id}
    mv ${id}_*.fq.gz ${id}

    echo "OK" > status.txt

    """

    stub:
    """
    mkdir ${id}
    touch ${id}/reads.fq.gz
    touch ${id}/report.txt
    touch ${id}/status.txt
    """

}