
process MERGEMPA {

    label 'process_single'

    input:
    val(sample_number)
    val(finish_number)
    path(abundance)
    path(samplesheet)

    output:
    path("mpa_report"),emit:"mpa_species"
    path("*.xls"),emit:"mpa_profile"

    when:
    sample_number == finish_number

    script:
    """

    merge_metaphlan_tables.py ${abundance} | metaphlan_fill_missing_species.pl - > abundance_table.txt

    col_reorder.pl abundance_table.txt 1 2 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_MetaPhlAn_abundance_table.xls

    parallel "mpa_tax_level.py -i ${params.pipeline_prefix}_MetaPhlAn_abundance_table.xls -l {1}; mv {1}.xls ${params.pipeline_prefix}_MetaPhlAn_{1}.xls" ::: phylum class order family genus species

    paste_mpa2sunburst.py ${params.pipeline_prefix}_MetaPhlAn_abundance_table.xls mpa.sunburst.txt

    # Delete result files before sorting.
    rm -f abundance_table.txt

    mkdir mpa_report
    cp ${params.pipeline_prefix}_MetaPhlAn_species.xls mpa_report/mpaspeciesT.xls
    cp ${params.pipeline_prefix}_MetaPhlAn_species.xls mpa_report/mpaspeciesTPCA.xls
    cp mpa.sunburst.txt mpa_report/mpaspecies.txt

    sed -i '1d' ${params.pipeline_prefix}_MetaPhlAn_abundance_table.xls
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS

    """

}
