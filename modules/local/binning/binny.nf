
process BINNY {
    
    tag "$id"
    
    label 'process_medium','error_ignore','process_long'

    input:
    tuple val(id),path(contigs),path(sorted_bam)

    output:
    tuple val(id),path("${id}/bins/",type:'dir'),emit:"bins"
    tuple val(id),path("binny.contigs2bin.tsv"),emit:"tsv"


    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.binny_options ?: ""
    
    """
    
    template_binny.py -f ${contigs} -b ${sorted_bam} -o \$PWD/${id} > ${id}.binny.yaml
    
    ${params.binny_path}/binny -l ${id}.binny.yaml -t ${task.cpus} ${options}

    finish=0
    if [ -d "${id}/bins" ] ; then      
        finish=\$((ls -1 "./${id}/bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then

        cd ${id}/bins
        i=1; for file in *; do new_name="binny_${id}_bin.\$i.fa"; mv "\$file" "\$new_name"; i=\$((i+1)); done
    
        cd ../..
        Fasta_to_Contig2Bin.sh -i ${id}/bins -e fa > binny.contigs2bin.tsv

        #sed -i 's%\\tbin%\\t${id}%g' binny.contigs2bin.tsv

    else
        touch binny.contigs2bin.tsv
    fi
    
    rm -rf ${id}/intermediary ${id}/contig_data.tsv.gz

    """

}
