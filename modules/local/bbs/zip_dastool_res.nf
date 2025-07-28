
process ZIPDASTOOLRES {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(all_eval),path(contig2bin)

    output:
    path("*tar.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    mkdir ${id}_DASTool
    cp ${all_eval} ${id}_DASTool
    cp ${contig2bin} ${id}_DASTool

    tar -zcf ${id}_DASTool.tar.gz ${id}_DASTool
    rm -rf ${id}_DASTool

    """
    stub:
    """
    touch ${id}_DASTool.tar.gz
    """

}
