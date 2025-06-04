
process RENAMECONTIGTAXO { 
    
    tag "$id"

    label 'process_low'
    
    input:
    tuple val(id),path(taxonomy),path(contig_map)

    output:
    tuple val(id), path("${id}_newID_contig_taxonomy.tsv"),emit:"newID_contig_taxonomy"
    

    when:
    task.ext.when == null || task.ext.when

    script:

    """

    rename_contig_taxonomy_contigID.py --map ${contig_map} --input ${taxonomy} --output ${id}_newID_contig_taxonomy.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$( python --version )
    END_VERSIONS


    """

    stub:
    """
    touch  ${id}_newID_contig_taxonomy.tsv

    """

}
