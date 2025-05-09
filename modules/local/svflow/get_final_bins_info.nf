
process GETBINSINFO {

    label 'process_low'

    input:
    path(bins_taxon)
    path(raw_bin_quality)
    path(all_bins_info)
    path(all_bins_rename_map)

    output:
    path("${params.pipeline_prefix}_final_bins_info.xls"),emit:"bins_info"

    when:
    task.ext.when == null || task.ext.when

    script:

    """

    SV_get_final_bins_info.py -t ${bins_taxon} -i ${raw_bin_quality} -b ${all_bins_info} -r ${all_bins_rename_map} -o ${params.pipeline_prefix}_final_bins_info.xls

    """
}
