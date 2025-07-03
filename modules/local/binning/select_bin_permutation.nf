process SELECTPERMUTATION {

    tag "$id"
    
    label 'process_single'

    input:
    tuple val(id), path(contigs2tsv), path(quality_reports), path(contig)

    output:
    tuple val(id), path("${id}_PBO_quality_report.tsv"), emit: "bin_qs", optional: true
    tuple val(id), path("${id}_PBO_best_bins", type: 'dir'), emit: "best_bin", optional: true
    tuple val(id), path("PermutationBest.contigs2bin.tsv"), emit: "contigs2bin", optional: true
    tuple val(id), path("${id}_PBO_BinsContigs.tsv"), emit: "BinsContigs", optional: true
    path("${id}_PBO_Tool_warning.log"), emit: "pbo_error_log", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Combine all quality reports into one
    ls ${quality_reports} | head -n 1 | xargs head -n 1 > all.checkm2.report
    ls ${quality_reports} | while read a; do
        tail -n +2 "\$a"
    done >> all.checkm2.report

    # Run permutation-based bin selection
    permutation_based_bin_seltection.py \\
        -t ${contigs2tsv} \\
        -c all.checkm2.report \\
        -f ${contig} \\
        -r ${id}_PBO_quality_report.tsv \\
        -b ${id}_PBO_best_bins \\
        -s PermutationBest.contigs2bin.tsv \\
        -m ${params.completeness} -n ${params.contamination}

    # Check if the expected output file was generated
    if [ -s PermutationBest.contigs2bin.tsv ]; then
        awk -F "\\t" '{print \$2"\\t"\$1"\\tPBO"}' PermutationBest.contigs2bin.tsv > ${id}_PBO_BinsContigs.tsv
    else
    cat <<- OUTLOG > ${id}_PBO_Tool_warning.log
    ========== Start at: \$(date) ==========
    ### Step: ${task.process}
    PBO task for sample ${id} failed.
    No bins were selected or generated. This may be due to all candidate bins failing the completeness/contamination thresholds.
    Please check input contig2bin file: ${contigs2tsv}, and quality reports in: ${quality_reports}
    ========== End at: \$(date) ==========
    OUTLOG
    fi
    """
}
