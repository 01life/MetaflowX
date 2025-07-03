
process SEMIBIN2 {
    
    tag "$id"

    label 'process_single','error_ignore'

    input:
    tuple val(id),path(contigs),path(sorted_bam)

    output:
    tuple val(id),path("${id}/output_bins/",type:'dir'), emit:"bins", optional: true
    tuple val(id),path("semibin2.contigs2bin.tsv"), emit:"tsv", optional: true
    tuple val(id),path("${id}_SemiBin2_BinsContigs.tsv"), emit:"BinsContigs", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.semibin2_options ?: ""
    
    """
    
    SemiBin2 single_easy_bin -p ${task.cpus} --input-fasta ${contigs} --input-bam ${sorted_bam} -o ${id} ${options} || echo "SEMIBIN2 task for sample ${id} failed ......" > ${id}.log

    finish=0
    if [ -d "${id}/output_bins" ] && [ ! -e "${id}.log" ]; then      
        finish=\$((ls -1 "./${id}/output_bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        
        cd ${id}/output_bins/
        #for file in *.fa; do mv "\$file" "${id}_\${file}"; done
        i=1; for file in *; do new_name="SemiBin2_${id}_bin.\$i.fa"; mv "\$file" "\$new_name"; i=\$((i+1)); done

        cd ../..

        Fasta_to_Contig2Bin.sh -i ${id}/output_bins -e fa > semibin2.contigs2bin.tsv

        awk -F "\\t" '{print\$2"\\t"\$1"\\tSemiBin2"}' semibin2.contigs2bin.tsv  > ${id}_SemiBin2_BinsContigs.tsv

        #sed -i 's%\\tSemiBin_%\\t${id}.%g' semibin2.contigs2bin.tsv
    
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin2: \$( SemiBin2 --version )
    END_VERSIONS

    """

}
