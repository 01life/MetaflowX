
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
    
    sed 1d ${samplesheet} | parallel -k --col-sep ',' "echo -e '{1}\\t{1}_abundance.xls'" > abundance.list
    
    merge_identical_rownames_tables.pl -table abundance.list -o ${params.pipeline_prefix}_merged_abundance.xls

    mkdir -p genesetAbundance_report
    ln -s ../${params.pipeline_prefix}_merged_abundance.xls genesetAbundance_report/genesetAbundance.xls


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | sed 's/perl //g')
    END_VERSIONS

    """
    stub:
    """
    mkdir -p genesetAbundance_report
    touch ${params.pipeline_prefix}_merged_abundance.xls
    """
}
