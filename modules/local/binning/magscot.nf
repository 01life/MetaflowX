
process MAGSCOT {
    
    tag "$id"
    
    label 'process_medium'

    input:
    tuple val(id),path(contigs),path(faa),path(bins_contig_csv)
    path(magscot_folder)

    output:
    tuple val(id),path("${id}_MAGScoT.scores.out"), emit: "scores", optional: true
    tuple val(id),path("${id}_MAGScoT.refined.contig_to_bin.out"), emit: "refined_contig_to_bin", optional: true
    tuple val(id),path("${id}_MAGScoT.refined.out"), emit: "refined", optional: true
    tuple val(id),path("${id}_MAGScoT_bins"), emit: "bins", optional: true
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def bin_list = bins instanceof List ? bins.join(",") : "$bins"
    def magscot_options = params.magscot_options ?: "" 
    """

    
    hmmsearch -o ${id}.hmm.tigr.out --tblout ${id}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${task.cpus} ${MAGScoT_folder}/hmm/gtdbtk_rel207_tigrfam.hmm ${faa}
    hmmsearch -o ${id}.hmm.pfam.out --tblout ${id}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${task.cpus} ${MAGScoT_folder}/hmm/gtdbtk_rel207_Pfam-A.hmm ${faa}

    cat ${id}.hmm.tigr.hit.out | grep -v "^#" | awk '{print \$1"\t"\$3"\t"\$5}' > ${id}.tigr
    cat ${id}.hmm.pfam.hit.out | grep -v "^#" | awk '{print \$1"\t"\$4"\t"\$5}' > ${id}.pfam
    cat ${id}.pfam ${id}.tigr > ${id}.hmm

    cat ${bins_contig_csv} > ${id}.contigs_to_bin_order.tsv
    Rscript ${MAGScoT_folder}/MAGScoT.R -i ${id}.contigs_to_bin_order.tsv --hmm ${id}.hmm -o ${id}_MAGScoT ${magscot_options}  || {
    
    cat <<-OUTLOG > ${id}_MAGScoT_error.txt

    ==========Start at : `date` ==========
    ### Step ${task.process}
    MAGScoT task for sample ${id} failed. It is spossible that no bins provided, please check contig2bin files: ${bin_list}. Alternatively, no bins found with scores greater than 0.5(default). This can occur if the input data is of poor quality or if the parameters used to run MAGScoT are not optimal.
    ==========End at : `date` ==========

    OUTLOG
    
    }

    #Process output during normal execution.

    split_contigs_to_bins.py -b ${id}_MAGScoT.refined.contig_to_bin.out -c ${contigs}   -o ${id}_MAGScoT_bins -p ${id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAGScoT: v1.1
    END_VERSIONS
    """

}
