
process TAXVAMB { 
    
    tag "$id"

    label 'process_medium'
    
    input:
    tuple val(id),path(contigs),path(depth),path(taxonomy)

    output:
    tuple val(id), path("${id}-taxvamb"),emit:"taxvamb"
    when:
    task.ext.when == null || task.ext.when

    script:
    def taxometer_options = params.taxometer_options ?: ""

    """
    vamb bin taxvamb --outdir ${id}-taxvamb --fasta ${contigs} --abundance_tsv ${depth} --taxonomy ${taxonomy} -p ${task.cpus}   ${taxometer_options}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Vamb: \$( vamb --version )
    END_VERSIONS

    """

    stub:
    """
    mkdir ${id}-taxvamb

    """

}
