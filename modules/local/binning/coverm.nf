process COVERM {

    tag "$method"

    label 'process_medium'

    input:
    val(sample_number)
    val(bam_number)
    tuple val(method),path(genomes)
    path(sorted_bam)
    path(samplesheet)

    output:
    path("${params.pipeline_prefix}_CoverM_bins_${method}.xls"),emit:"abundance"
    path("coverm_${method}_report"),emit:"report"

    when:
    sample_number == bam_number && method != 'coverage_histogram'

    script:
    def options = params.coverm_options ?: ""
    """
    coverm genome -b ${sorted_bam} -d ${genomes} -m ${method} -t ${task.cpus} -o ${params.pipeline_prefix}_CoverM_bins_${method}_raw.xls ${options}

    standardize_coverm_outfile.py -i ${params.pipeline_prefix}_CoverM_bins_${method}_raw.xls -o ${params.pipeline_prefix}_CoverM_bins_${method}_rename.xls -m ${method}

    col_reorder.pl ${params.pipeline_prefix}_CoverM_bins_${method}_rename.xls 1 1 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_CoverM_bins_${method}.xls

    mkdir coverm_${method}_report
    cp ${params.pipeline_prefix}_CoverM_bins_${method}_raw.xls coverm_${method}_report/bin_${method}.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/coverm //g')
    END_VERSIONS

    """

}
