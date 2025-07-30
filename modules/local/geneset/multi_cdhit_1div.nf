
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
    #mkdir div_tmp
    psort2div.pl --input ${allcds} --output div_tmp --cpu ${task.cpus} -prefix all.cds.fa.div -number ${split_num}

    #ln -s ${allcds} all.cds.fa
    #cd-hit-div -i all.cds.fa -o div_tmp/all.cds.fa.div -div ${split_num}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit-est 2>&1 | grep 'CD-HIT version' | sed 's/CD-HIT //')
    END_VERSIONS

    """

    stub:
    """
    mkdir -p div_tmp
    touch div_tmp/all.cds.fa.div-1
    touch div_tmp/all.cds.fa.div-2
    touch div_tmp/all.cds.fa.div-3
    touch div_tmp/all.cds.fa.div-4
    touch div_tmp/all.cds.fa.div-5
    touch div_tmp/all.cds.fa.div-6
    touch div_tmp/all.cds.fa.div-7
    touch div_tmp/all.cds.fa.div-8
    touch div_tmp/all.cds.fa.div-9
    touch div_tmp/all.cds.fa.div-10
    touch div_tmp/all.cds.fa.div-11
    """

}
