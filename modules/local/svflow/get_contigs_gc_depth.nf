
process GETGCANDDEPTH {

    label 'process_low'

    input:
    path(contig_info)
    path(depth_list)
    path(bin_list)

    output:
    path("${params.pipeline_prefix}_contigs_gc_depth.xls"),emit:"info"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    for f in *_depth.xls; do echo -e "\${f%_depth.xls}\\t\$f"; done > all.depth.list
    ls folder*/*fa |awk -F "/" '{print\$NF"\\t"\$0}' |sed 's/.fa\\t/\\t/g' > bin.fa.list

    SV_get_contig_info.py -i ${contig_info} -d all.depth.list -b bin.fa.list -o ${params.pipeline_prefix}_contigs_gc_depth.xls


    """

}
