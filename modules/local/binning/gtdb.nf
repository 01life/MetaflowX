
process GTDB {

    tag "$baseName"

    label 'process_high'

    input:
    tuple val(baseName),path(genomes)
    path(gtdbtk_db)
    path(mash_db)
    path(gtdb_archaeal_metadata)
    path(gtdb_bacterial_metadata)

    output:
    path("gtdb_output_${baseName}")
    path("${baseName}*gtdbtk.tsv"),emit:"gtdbtk_summary"
    path("gtdb2ncbi_${baseName}.txt"),emit:"gtdb2ncbi"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.gtdbtk_classify_wf_options ?: ""
    """

    export GTDBTK_DATA_PATH="${gtdbtk_db}"
    
    gtdbtk classify_wf --genome_dir ${genomes} --out_dir gtdb_output_${baseName} --cpus ${task.cpus} --mash_db ${mash_db} ${options}

    gtdb_to_ncbi_majority_vote.py \
        --ar53_metadata_file ${gtdb_archaeal_metadata} \
        --bac120_metadata_file ${gtdb_bacterial_metadata} \
        --gtdbtk_output_dir gtdb_output_${baseName} \
        --output_file gtdb2ncbi_${baseName}.txt

    if [ -e "gtdb_output_${baseName}/gtdbtk.bac120.summary.tsv" ]; then
        if [ -e "gtdb_output_${baseName}/gtdbtk.ar53.summary.tsv" ]; then
            #bac120 and ar53

            # Get the header row.
            find "gtdb_output_${baseName}" -maxdepth 1 -name "gtdbtk.*.summary.tsv" | head -n 1 | xargs head -n 1 > "${baseName}.bac120.ar53.gtdbtk.tsv"

            # Merge all non-header rows.
            find "gtdb_output_${baseName}" -maxdepth 1 -name "gtdbtk.*.summary.tsv" | while read a; do
                tail -n +2 "\$a"
            done >> "${baseName}.bac120.ar53.gtdbtk.tsv"

        #only bac120
        else
            cp gtdb_output_${baseName}/gtdbtk.bac120.summary.tsv ${baseName}.bac120.gtdbtk.tsv
        fi
    else
        #only ar53
        if [ -e "gtdb_output_${baseName}/gtdbtk.ar53.summary.tsv" ]; then
            cp gtdb_output_${baseName}/gtdbtk.ar53.summary.tsv ${baseName}.ar53.gtdbtk.tsv
        #empty result
        else
            touch ${baseName}.empty.gtdbtk.tsv
        fi
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk --version | sed -n 1p | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS

    """
}
