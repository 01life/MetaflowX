
process SEMIBIN2 {
    
    tag "$id"

    label 'process_single','error_ignore'

    input:
    tuple val(id),path(contigs),path(sorted_bam)

    output:
    tuple val(id),path("${id}/output_bins/",type:'dir'),emit:"bins"
    tuple val(id),path("semibin2.contigs2bin.tsv"),emit:"tsv"


    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.semibin2_options ?: ""
    
    """
    
    SemiBin2 single_easy_bin \
        -p ${task.cpus} \
        --input-fasta ${contigs} \
        --input-bam ${sorted_bam} \
        -o ${id} ${options}

    finish=0
    if [ -d "${id}/output_bins" ] ; then      
        finish=\$((ls -1 "./${id}/output_bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        
        cd ${id}/output_bins/
        for file in *.fa; do mv "\$file" "${id}_\${file}"; done
        cd ../..

        Fasta_to_Contig2Bin.sh -i ${id}/output_bins -e fa > semibin2.contigs2bin.tsv

        #sed -i 's%\\tSemiBin_%\\t${id}.%g' semibin2.contigs2bin.tsv
    
    else
        touch semibin2.contigs2bin.tsv
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin2: \$( SemiBin2 --version )
    END_VERSIONS

    """

}
