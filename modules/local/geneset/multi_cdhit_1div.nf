
process CDHITDIV {

    label 'process_single'

    input:
    path(allcds)
    path(task_num)

    output:
    path("div_tmp/all.cds.fa.div-*"),emit:"div"
    path("div_tmp"),emit:"div_tmp"

    script:
    def split_num = task_num.getSimpleName().toInteger()
    
    """
    mkdir div_tmp
    ln -s ${allcds} all.cds.fa
    cd-hit-div -i all.cds.fa -o div_tmp/all.cds.fa.div -div ${split_num}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit-est 2>&1 | grep 'CD-HIT version' | sed 's/CD-HIT //')
    END_VERSIONS

    """

}
