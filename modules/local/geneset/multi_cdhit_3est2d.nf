
process CDHITEST2D {

    tag "$div_num"

    label 'process_single'

    input:
    tuple val(div_num),path(div_tmp),path(div0)
    
    output:
    path("clstr_tmp/*.clstr"),emit:"clstr"
    path("clstr_tmp/all.cds.fa.div-*-o"),emit:"div_o"
    
    script:
    def cdhitoptions = params.cdhit_options ?: ""

    
    """
    mkdir est_tmp

    ln -s \$(readlink -f ${div0}) ./est_tmp

    for file in ${div_tmp}/*; do
    ln -s \$(readlink -f \$file) ./est_tmp/\$(basename \$file)
    done


    i=${div_num}
    indiv=est_tmp/all.cds.fa.div

    mkdir clstr_tmp

    # Create the file for segment \$i
    idb="\${indiv}-\${i}"
    idblog="clstr_tmp/all.cds.fa.div-\${i}.log"

    # Compare to previous segments
    for (( j=0; j<i; j++ ))
    do
        jdb="\${indiv}-\${j}-o"
        idbo="clstr_tmp/all.cds.fa.div-\${i}.vs.\${j}"
        cd-hit-est-2d -i \$jdb -i2 \$idb -o \$idbo -T ${task.cpus} ${cdhitoptions} >> \$idblog

        idb=\$idbo
    done

    # Self-comparing
    cd-hit-est -i \$idb -o clstr_tmp/all.cds.fa.div-\${i}-o -T ${task.cpus} ${cdhitoptions} >> \$idblog


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit-est 2>&1 | grep 'CD-HIT version' | sed 's/CD-HIT //')
    END_VERSIONS

    """

}
