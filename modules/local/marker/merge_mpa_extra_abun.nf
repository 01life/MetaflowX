
process MERGEMPAEXTRAABUN {

    tag "$method"

    label 'process_single'

    input:
    tuple val(method),path(abundance)

    output:
    path("${params.pipeline_prefix}_MetaPhlAn_${method}.xls"),emit:"extra_profile"

    when:
    method == "rel_ab_w_read_stats"

    script:
    """

    merge_metaphlan_tables.py ${abundance} > ${params.pipeline_prefix}_MetaPhlAn_${method}.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS

    """

    stub:
    """
    touch ${params.pipeline_prefix}_MetaPhlAn_${method}.xls
    """

}
