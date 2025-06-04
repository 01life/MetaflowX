
process ORDERCOMMONID { 
    
    tag "$id"

    label 'process_low'
    
    input:
    tuple val(id),path(contigs),path(depth),path(taxonomy)

    output:
    tuple val(id), path("${id}_contigs_common.fa"),path("${id}_depth_common.txt"),path("${id}_contig_taxonomy_common.tsv"),emit:"commoninput"
    

    when:
    task.ext.when == null || task.ext.when

    script:

    """

    awk -v s="${id}" 'BEGIN{OFS="\\t"} NR==1{print "contigname", s; next} {print \$1, \$3}' ${depth} > abundances.tsv

    extract_common_ids4taxometer.py abundances.tsv ${taxonomy} ${contigs}  --out1 ${id}_depth_common.txt --out2 ${id}_contig_taxonomy_common.tsv --out3 ${id}_contigs_common.fa


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$( python --version )
    END_VERSIONS


    """

    stub:
    """
    touch  ${id}_depth_common.txt
    touch  ${id}_contig_taxonomy_common.tsv 
    touch  ${id}_contigs_common.fa

    """

}
