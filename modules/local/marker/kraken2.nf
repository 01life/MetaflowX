
process KRAKEN2 {
    
    tag "$id"

    label 'process_high'
    
    input:
    tuple val(id),path(reads)
    path(kraken2_db)

    output:
    path("${id}*.xls")
    path("${id}_bracken_*_mpa.xls"),emit:"mapping"
    path("${id}_bracken_species_mpa.xls"),emit:"species"

    when:
    task.ext.when == null || task.ext.when

    script:
    def paired = params.single_end ? "" : "--paired"
    def kraken_options = params.kraken2_options ?: ""
    def bracken_options = params.bracken_options ?: ""
    """

    kraken2 --db ${kraken2_db} --threads ${task.cpus} --report ${id}_kreport.xls ${paired} ${reads} ${kraken_options} > ${id}_kraken.xls

    parallel --xapply "bracken -d ${kraken2_db} -i ${id}_kreport.xls -o ${id}_bracken_{2}.xls -l {1} -t ${task.cpus} ${bracken_options}; kreport2mpa.py --no-intermediate-ranks --display-header -r ${id}_kreport_bracken_{2}.xls -o ${id}_bracken_{2}_mpa.xls " ::: D P C O F G S ::: domains phylums classes orders families genuses species

    rm -rf ${id}_kraken.xls
        
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
