
process GETBINSTAXON {

    label 'process_low'

    input:
    path(gtdb_summary)
    path(bins_st_rel_abun)
    
    output:
    path("${params.pipeline_prefix}_bins_taxonomic.xls"),emit:"taxon"
    path("${params.pipeline_prefix}*.xls"),emit:"profile"

    when:
    gtdb_summary.size() > 0

    script:

    """
        
    SV_gettaxon.py ${gtdb_summary} ${params.pipeline_prefix}_bins_taxonomic.xls

    parallel "SV_bin_tax_level.py -i ${params.pipeline_prefix}_bins_taxonomic.xls -a ${bins_st_rel_abun} -l {1} -p ${params.pipeline_prefix} -o ./" ::: phylum class order family genus species

    SV_get_lefse_abundance_table.py -i ${gtdb_summary} -a ${bins_st_rel_abun} -o ${params.pipeline_prefix}_bins_abundance_table.xls

    SV_get_mergesankey.py ${bins_st_rel_abun} ${params.pipeline_prefix}_bins_taxonomic.xls ${params.pipeline_prefix}_sankey.xls

    """

}
