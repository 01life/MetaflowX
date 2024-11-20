process GETBINMASHFQ {

    label 'process_low'

    input:
    val(outdir)
    path(clean_reads)
    path(filter_bin_info)
    path(bin_taxon)
    
    output:
    path("${params.pipeline_prefix}_bin_mash_fq.txt"),emit:"bin_mash_fq"

    when:
    task.ext.when == null || task.ext.when

    script:
    def remove_samples = params.remove_samples ? "--remove_samples ${params.remove_samples}" : ""
    def options = "-q ${params.min_quality_score} -a ${params.bin_min_abundance} -p ${params.bin_min_popularity} ${remove_samples}"

    """

    get_input_bin_reassembly_V2.py -b ${outdir} \\
        --fq_paths_file ${clean_reads} \\
        --gtdb_genome_paths_file ${params.genome_paths} \\
        --bin_quality_file  ${filter_bin_info} \\
        --prefix ${params.pipeline_prefix} \\
        -o ${params.pipeline_prefix}_bin_mash_fq.txt ${options}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS


    """

}
