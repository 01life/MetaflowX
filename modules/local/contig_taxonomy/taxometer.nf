
process TAXOMETER { 
    
    tag "$id"

    label 'process_medium'
    
    input:
    tuple val(id),path(contigs),path(depth),path(taxonomy)

    output:
    tuple val(id), path("${id}-taxometer"),emit:"taxometer"
    // tuple val(id), path("${id}_taxometer_0.95"),emit:"taxometer95"
    // tuple val(id), path("${id}_taxometer_0.95/${id}_species.tsv"),emit:"taxometer95Species"
    when:
    task.ext.when == null || task.ext.when

    script:
    def kraken_options = params.kraken2_contig_anno_options ?: ""
    """
    cut -f 1,2  ${taxonomy} >  ${id}.taxonomy.tmp
    vamb taxometer --outdir ${id}-taxometer --fasta ${contigs} --abundance_tsv ${depth} --taxonomy ${id}.taxonomy.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Vamb: \$( vamb --version )
    END_VERSIONS

    """

    stub:
    """
    mkdir ${id}-taxometer

    """

}
