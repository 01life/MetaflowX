process GETBINSAMPLE {

    label 'process_low'

    input:
    path(bin_summary)
    path(bin_count_abundance)
    path(bin_mean_abundance)

    output:
    path("${params.pipeline_prefix}_pick4optimize_bin_sample.txt") ,emit:"binsample"

    when:
    task.ext.when == null || task.ext.when

    script:
    def get_bin_assembly_options = params.get_bin_assembly_options ?: ""

   
    """
    bra_get_reassembly_bin_sample.py \\
        -i ${bin_summary} \\
        -c ${bin_count_abundance} \\
        -m ${bin_mean_abundance} \\
        --prefix ${params.pipeline_prefix} \\
        --gtdb_genome_paths_file ${params.genome_paths}  \\
        ${get_bin_assembly_options}
        

    """
}
