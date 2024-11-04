
process ESTDIV {

    label 'process_single'

    input:
    path(div)

    output:
    path("${div}-o"),emit:"div0"

    script:
    def options = params.cdhit_options ?: ""

    """
    
    cd-hit-est -i ${div} -o  ${div}-o  -T  ${task.cpus}  ${options}   >>  ${div}-o.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit-est 2>&1 | grep 'CD-HIT version' | sed 's/CD-HIT //')
    END_VERSIONS

    """

}
