
process MERGEBINABUNTAXON {

    label 'process_single'

    input:
    path(gtdb_summary)
    path(bins_rel_abun)
    
    output:
    path("${params.pipeline_prefix}*.csv"),emit:"taxon"

    when:
    gtdb_summary.size() > 0

    script:

    """

    merge_bin_abundace_taxonomy.py -t ${gtdb_summary} -p ${bins_rel_abun} -q ${params.pipeline_prefix}

    """
    stub:
    """
    touch ${params.pipeline_prefix}_bins_taxonomy.csv
    touch ${params.pipeline_prefix}_bins_rel_abundance.csv
    """
}
