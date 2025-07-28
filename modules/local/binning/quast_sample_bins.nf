process QUAST_SAMPLE_BINS {

    tag "$id"
    
    label 'process_medium'

    input:
    tuple val(id), path(bins_folder)

    output:
    path "${id}.QUAST.tar.gz", emit: quast_folder
    tuple val(id), path("${id}_Bins_merged_qsast_summary.tsv"), emit: quast_bin_summaries
    // path "versions.yml", emit: versions

    script:
    def metaquastOptions = params.metaquast_options ?: ""

    """
    mkdir -p QUAST

    for contig in ${bins_folder}/*.fa; do
        fname=\$(basename "\$contig" .fa)
        outdir="QUAST/\$fname"
        mkdir -p "\$outdir"
        metaquast.py --threads ${task.cpus} -l "\$fname" "\$contig" -o "\$outdir" ${metaquastOptions}
    done

    awk 'FNR==1 && NR!=1 {next;}{print}'  QUAST/*/transposed_report.tsv > ${id}_merged_qsast_summary.tmp 
        
    merge_quast_sample_bins.py  --summaries ${id}_merged_qsast_summary.tmp  --result_dir ./QUAST  --output_file ${id}_Bins_merged_qsast_summary.tsv

    tar -zcf ${id}.QUAST.tar.gz QUAST

    rm -rf ${id}_merged_qsast_summary.tmp QUAST


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        metaquast: \$(metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p QUAST
    touch ${id}.QUAST.tar.gz
    touch ${id}_Bins_merged_qsast_summary.tsv
    """
}
