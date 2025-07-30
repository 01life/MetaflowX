
process VAMBBIN {
    
    tag "$id"
    
    label 'process_single','error_ignore'
    
    input:
    tuple val(id),path(contigs),path(bin_depth)

    output:
    tuple val(id),path("${id}/VambBins",type:'dir'), emit:"bins", optional: true
    tuple val(id),path("vamb.contigs2bin.tsv"), emit:"tsv", optional: true
    tuple val(id),path("${id}_Vamb_BinsContigs.tsv"), emit:"BinsContigs", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.vamb_options ?: ""
    """


    cut -f 1,3 ${bin_depth}  | awk -v sname="${id}" 'NR==1 {print "contigname\\t"sname; next} {print}' > abundances.tsv

    vamb bin default --outdir ${id} --fasta ${contigs} --abundance_tsv abundances.tsv -p  ${task.cpus}  ${options} || echo "vamb task for sample ${id} failed ......" > ${id}.log


    finish=0
    if [ -d "${id}" ] && [ ! -e "${id}.log" ]; then      
        finish=\$((ls -1 "${id}/bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        cd ${id}/bins
            i=1; for file in *; do new_name="vamb_${id}_bin.\$i.fa"; mv "\$file" "\$new_name"; i=\$((i+1)); done

        cd ..
            Fasta_to_Contig2Bin.sh -i bins -e fa > vamb.contigs2bin.tsv

            awk -F "\\t" '{print\$2"\\t"\$1"\\tVamb"}' vamb.contigs2bin.tsv  > ${id}_Vamb_BinsContigs.tsv

        mv bins VambBins
        mv vamb.contigs2bin.tsv ../
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Vamb: \$( vamb --version )
    END_VERSIONS

    """

    stub:
    """
    mkdir -p ${id}/VambBins
    touch ${id}/VambBins/vamb_${id}_bin.1.fa
    touch vamb.contigs2bin.tsv
    touch ${id}_Vamb_BinsContigs.tsv
    """
}
