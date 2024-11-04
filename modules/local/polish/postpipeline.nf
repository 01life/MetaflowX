process POSTPIPELINE {
    
    tag "$id"
    
    label 'process_low'

    input:
    path(report)
    tuple val(id),path(reads)

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

}