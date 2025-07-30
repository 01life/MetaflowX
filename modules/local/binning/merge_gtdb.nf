
process MERGEGTDB {

    label 'process_single'

    input:
    val(summary_linecount)
    path(gtdbtk_summary)
    path(gtdb2ncbi)
    path(qs_quality_report)
    path(bins_rename_map)
    path(sample_bin_map)


    output:
    path("gtdbtk.summary.tsv"), emit: "gtdb_summary"
    path("gtdb_report"), emit: "gtdb_report", optional: true
    path("bin_QS_taxonomy_summary.xls"), emit: "bin_QS_taxonomy", optional: true

    when:
    summary_linecount > 1

    script:
    """
    cp ${gtdbtk_summary} gtdbtk.summary.tsv
    cp ${gtdb2ncbi} gtdbtk.taxonomy2ncbi.summary.tsv

    paste_gtbk2treemap.py gtdbtk.summary.tsv gtdb2treemap.txt

    mkdir gtdb_report
    cp gtdb2treemap.txt gtdb_report/gtdb.txt

    merge_bin_info.py -g ${gtdb2ncbi} \
        -c ${qs_quality_report}  \
        -m ${bins_rename_map}  \
        -s ${gtdbtk_summary} \
        -d ${sample_bin_map} \
        -o bin_QS_taxonomy_summary.xls

    cp bin_QS_taxonomy_summary.xls gtdb_report/bintable.xls

    """
    
    stub:
    """
    touch gtdbtk.summary.tsv
    touch gtdb2treemap.txt
    mkdir gtdb_report
    touch gtdb_report/gtdb.txt
    touch bin_QS_taxonomy_summary.xls
    touch gtdb_report/bintable.xls
    """

}
