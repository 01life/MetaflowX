
process TAXOMETERFILTER { 
    
    tag "$id"

    label 'process_low'
    
    input:
    tuple val(id),path(taxometer)

    output:
    tuple val(id), path("${id}_taxometer_0.95"),emit:"taxometer95"
    tuple val(id), path("${id}_taxometer_0.95/${id}_species.tsv"),emit:"taxometer95Species"
    tuple val(id), path("${id}_taxometer_0.95/${id}_deepest_level.tsv"),emit:"taxometer95deepest"

    when:
    task.ext.when == null || task.ext.when

    script:
    def kraken_options = params.kraken2_contig_anno_options ?: ""
    """

    filter_taxometer_score.py -i  ${taxometer} -o ${id}_taxometer_0.95 -p ${id} --score_threshold 0.95

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$( python --version )
    END_VERSIONS

    """

    stub:
    """
    mkdir ${id}-taxometer

    """

}
