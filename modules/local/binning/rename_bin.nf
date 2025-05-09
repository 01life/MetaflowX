
process RENAMEBIN {

    label 'process_low'

    input:
    path(final_bins)

    output:
    path("HQUniqueBins/folder*"),emit:"folder"
    path("HQUniqueBins/folder*/*.fa"),emit:"genomes"
    path("binInfo_report"),emit:"binInfo"
    path("${params.pipeline_prefix}_HQ_unique_bins_rename_map.xls"),emit:"name_map"
    path("${params.pipeline_prefix}_HQ_unique_bins_info.xls"),emit:"bins_info"
    path("*.txt"),emit:"count"
    path("HQUniqueBins/folder*"),emit:"list"

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    
    bin_rename.py -i ./ -o HQUniqueBins -p ${params.pipeline_prefix}

    split_folder.py HQUniqueBins ${params.gtdb_bin_chunk_size} 
    mv HQUniqueBins/*.txt .

    ls \$PWD/HQUniqueBins/folder*/*fa |awk -F "/" '{print\$NF"\\t"\$0}' |sed 's/.fa\\t/\\t/g' > bin.fa.list

    bin_stat.py bin.fa.list ${params.pipeline_prefix}_HQ_unique_bins_info.xls

    mkdir binInfo_report
    cp ${params.pipeline_prefix}_HQ_unique_bins_info.xls binInfo_report/binInfo.xls


    """
}
