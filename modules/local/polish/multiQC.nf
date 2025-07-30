
process MULTIQC {

    label 'process_low'

    input:
    val(type)
    path(multiqc_files)

    output:
    path("${type}.html"), emit:"html"
    path("${type}_data/"),emit:"data"

    script:

    """
    
    multiqc . -n ${type}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    
    """

    stub:
    """
    touch ${type}.html
    mkdir ${type}_data/
    """

}
