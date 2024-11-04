
process MERGEGENEPROFILE {

    label 'process_single'

    input:
    val(sample_number)
    val(finish_number)
    path(abundance)
    path(samplesheet)

    output:
    path("${params.pipeline_prefix}_*.xls"),emit:"merge"
    path("genesetAbundance_report"),emit:"geneset_abundance_report"

    when:
    sample_number == finish_number

    script:
    """
    
    echo -e "id\\tfile" > abundance.list
    ls -d \$PWD/*_abundance.xls | awk -v pwd="\$PWD" '{sub(/.*\\//, "", \$0); sub(/_abundance\\.xls\$/, "", \$0); print \$0, pwd "/" \$0 "_abundance.xls"}' OFS='\\t' >> abundance.list

    merge.py -i abundance.list -p raw

    col_reorder.pl raw*.xls 1 1 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_merged_abundance.xls

    #Delete result files before sorting.
    rm -f raw*.xls

    mkdir genesetAbundance_report
    cp ${params.pipeline_prefix}_*.xls genesetAbundance_report/genesetAbundance.xls 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS

    """

}
