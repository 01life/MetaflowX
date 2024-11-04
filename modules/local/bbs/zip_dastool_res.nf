
process ZIPDASTOOLRES {
    
    label 'process_low'

    input:
    path(all_eval)
    path(contig2bin)
    path(contig_list)

    output:
    path("DASTool.tar.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    mkdir DASTool
    cp ${all_eval} DASTool/
    cp ${contig2bin} DASTool/

    cd DASTool
    awk -F "\\t" '{print "mkdir "\$1"\\n mv "\$1"*.t* ./"\$1}' ../${contig_list} |sh
    cd ..

    tar -zcf DASTool.tar.gz DASTool
    rm -rf DASTool

    """

}
