
process KRKEN2CONTIGTAXO  {
    
    tag "$id"

    label 'process_high'
    
    input:
    tuple val(id),path(contigs)
    path(kraken2_db)

    output:
    tuple val(id), path("${id}_kreport.xls"),emit:"kraken_report"
    tuple val(id), path("${id}_output.txt"),emit:"kraken_output"
    tuple val(id), path("${id}_contig_Kraken2_taxonomy.txt"),emit:"kraken_taxonomy"

    when:
    task.ext.when == null || task.ext.when

    script:
    def kraken_options = params.kraken2_contig_anno_options ?: ""
    """

    kraken2 --db ${kraken2_db} --use-names --memory-mapping  --threads ${task.cpus}  --report ${id}_kreport.xls --output ${id}_output.txt  ${contigs} ${kraken_options} 

    parallel --xapply "bracken -d ${kraken2_db} -i ${id}_kreport.xls -o ${id}_bracken_{2}.xls -l {1} ${bracken_options}; kreport2mpa.py --no-intermediate-ranks --display-header -r ${id}_kreport_bracken_{2}.xls -o ${id}_bracken_{2}_mpa.xls " ::: D P C O F G S ::: domains phylums classes orders families genuses species

    translate_Kraken2_annotate_taxonomy.py --taxonomy ${id}_bracken_species_mpa.xls --kraken ${id}_output.txt  --output ${id}_contig_Kraken2_taxonomy.tsv


        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS

    """

    stub:
    """
    touch ${id}_bracken_species_mpa.xls

    """

}
