
process SBG2GTDB {

    tag "$id"

    label 'process_single'

    input:
    tuple val(id),path(profile)
    path(mpa_db)

    output:
    path("${id}_sgb2gtdb.xls"),emit:"profile"

    when:
    params.mpa_index in ["mpa_vJan21_CHOCOPhlAnSGB_202103"]

    script:
    """

    sgb_to_gtdb_profile.py -i ${profile} -d ${mpa_db}/${params.sgb2gtdb_index} -o ${id}_sgb2gtdb.xls 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS

    """

    stub:
    """
    touch ${id}_sgb2gtdb.xls
    """
}
