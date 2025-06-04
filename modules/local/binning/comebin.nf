
process COMEBIN {
    
    tag "$id"
    
    label 'process_medium','error_ignore','process_long'

    input:
    tuple val(id),path(contigs),path(sorted_bam)

    output:
    tuple val(id),path("${id}/comebin_res/comebin_res_bins/",type:'dir'), emit:"bins", optional: true
    tuple val(id),path("comebin.contigs2bin.tsv"), emit:"tsv", optional: true
    tuple val(id),path("${id}_COMEBin_BinsContigs.tsv"), emit:"BinsContigs", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.comebin_options ?: ""

    """
    
    mkdir bam_file
    mv ${sorted_bam} bam_file

    run_comebin.sh -a ${contigs} -p bam_file -t ${task.cpus} -o ${id} ${options} || echo "COMEBIN task for sample ${id} failed ......" > ${id}.log
    
    finish=0
    if [ -d "${id}/comebin_res/comebin_res_bins" ] && [ ! -e "${id}.log" ] ; then      
        finish=\$((ls -1 "./${id}/comebin_res/comebin_res_bins") | wc -l)
        
    fi

    if [ \$finish -gt 0 ]; then
        cd ${id}/comebin_res/comebin_res_bins
        i=1; for file in *; do new_name="comebin_${id}_bin.\$i.fa"; mv "\$file" "\$new_name"; i=\$((i+1)); done
    
        cd ../../..

        Fasta_to_Contig2Bin.sh -i ${id}/comebin_res/comebin_res_bins/ -e fa > comebin.contigs2bin.tsv

        awk -F "\\t" '{print\$2"\\t"\$1"\\tCOMEBin"}' comebin.contigs2bin.tsv > ${id}_COMEBin_BinsContigs.tsv

        #sed -i 's%\\tbin%\\t${id}%g' comebin.contigs2bin.tsv
    
    fi

    rm -rf ${id}/data_augmentation ${id}/comebin_res/cluster_res ${id}/comebin_res/checkpoint_0200.pth.tar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh 2>&1 | sed -n '2p' | awk -F':' '{print \$2}')
    END_VERSIONS 

    """

}
