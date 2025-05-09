
process BUILDTREE {

    tag "$baseName"

    label 'process_medium'

    input:
    path(count)
    tuple val(baseName),path(genomes)

    output:
    path("gtdb_infer/infer/intermediate_results/gtdbtk.unrooted.tree"),emit:"tree"

    when:
    count.getSimpleName().toInteger() <= params.gtdb_bin_chunk_size

    script:
    def options1 = params.gtdbtk_identify_options ?: ""
    def options2 = params.gtdbtk_align_options ?: ""
    def options3 = params.gtdbtk_infer_options ?: ""

    """
    export GTDBTK_DATA_PATH="${params.gtdbtk_db}"
    
    gtdbtk identify --genome_dir ${genomes} --out_dir gtdb_identify --cpus ${task.cpus} ${options1}

    gtdbtk align --identify_dir gtdb_identify --out_dir gtdb_align --cpus ${task.cpus} ${options2}

    gtdbtk infer --msa_file gtdb_align/align/gtdbtk.bac120.user_msa.fasta.gz --out_dir gtdb_infer --cpus ${task.cpus} ${options3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version | sed -n 1p | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS

    """
}
