
process BRASEMIBIN2 {
    
    tag "$id-$sample"

    label 'process_single','error_ignore'

    input:
    tuple val(id),val(sample),path(contigs),path(sorted_bam)

    output:
    // tuple val(id),path("${id}/output_recluster_bins/*.fa"),emit:"fa"
    tuple val(id),val(sample),path("${id}/output_bins/",type:'dir'),emit:"bins"
    tuple val(id),val(sample),path("${id}_${sample}_semibin2.contigs2bin.tsv"),emit:"tsv"


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

    Fasta_to_Contig2Bin.sh -i ${id}/output_bins -e fa > ${id}_${sample}_semibin2.contigs2bin.tsv

    sed -i 's%\\tSemiBin_%\\t${id}.%g' ${id}_${sample}_semibin2.contigs2bin.tsv

    #touch ${id}_${sample}_semibin2.contigs2bin.tsv
    #mkdir -p ${id}/output_bins/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin2: \$(echo \$(SemiBin --version 2>&1) | sed 's/SemiBin //g' )
    END_VERSIONS

    """

}
