
process MERGEKRAKEN2 {

    label 'process_single'

    input:
    val(sample_number)
    val(finish_number)
    path(mapping)
    path(samplesheet)

    output:
    path("*.xls")
    path("kraken_report"),emit:"kraken2_species"

    when:
    sample_number == finish_number

    script:
    """
    
    ls *_bracken_*_mpa.xls |awk -F "_bracken_" '{print\$1"\t"\$0}' > tmp_mpa.txt

    parallel --xapply "grep {1} tmp_mpa.txt > {1}.list.mpa; merge_multi_matrix_XY_0.py {1}.list.mpa ${params.pipeline_prefix}_Kraken2_{2}_raw.txt; col_reorder.pl ${params.pipeline_prefix}_Kraken2_{2}_raw.txt 1 1 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_Kraken2_{2}.xls" ::: domains phylums classes orders families genuses species ::: domain phylum class order family genus species

    # Delete result files before sorting.
    rm -rf ${params.pipeline_prefix}_Kraken2_*_raw.txt

    paste_mpa2sunburst.py ${params.pipeline_prefix}_Kraken2_species.xls kraken.sunburst.txt
    
    mkdir kraken_report
    cp ${params.pipeline_prefix}_Kraken2_species.xls kraken_report/krakenspeciesT.xls
    cp ${params.pipeline_prefix}_Kraken2_species.xls kraken_report/krakenspeciesTPCA.xls
    cp kraken.sunburst.txt kraken_report/krakenspecies.txt
    
    """
    stub:
    """
    touch ${params.pipeline_prefix}_Kraken2_species.xls
    touch kraken.sunburst.txt
    mkdir kraken_report

    """ 

}
