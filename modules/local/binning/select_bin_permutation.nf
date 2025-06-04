process  SELECTPERMUTATION {
    
    tag "$id"
    
    label 'process_single'

    input:
    tuple val(id),path(contigs2tsv),path(quality_reports),path(contig)

    output:
    tuple val(id),path("${id}_PBO_quality_report.tsv"),emit:"bin_qs"
    tuple val(id),path("${id}_PBO_best_bins",type:'dir'),emit:"best_bin"
    tuple val(id),path("PermutationBest.contigs2bin.tsv"),emit:"contigs2bin"
    tuple val(id),path("${id}_PBO_BinsContigs.tsv"), emit:"BinsContigs", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    
    """

    ls ${quality_reports} | head -n 1 | xargs head -n 1 > all.checkm2.report

    ls ${quality_reports}| while read a; do
                tail -n +2 "\$a"
            done >> all.checkm2.report

    permutation_based_bin_seltection.py  -t ${contigs2tsv} -c all.checkm2.report -f ${contig} -r ${id}_PBO_quality_report.tsv -b ${id}_PBO_best_bins -s PermutationBest.contigs2bin.tsv -m ${params.completeness} -n ${params.contamination}

    awk -F "\\t" '{print\$2"\\t"\$1"\\tPBO"}'  PermutationBest.contigs2bin.tsv  > ${id}_PBO_BinsContigs.tsv

    
    """

}