process FASTQSTAT {
    
    tag "$id"

    label 'process_low'

    input:
    tuple val(id),path(raw_info),path(rmhost_log),path(clean_reads)

    output:
    path("${id}_reads_stat.txt"),emit:"reads_stat"
    path("${id}_clean_reads.txt"),emit:"clean_reads"

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_path = params.single_end ? "\$(readlink -f ${id}_clean.fq.gz)" : "\$(readlink -f ${id}_clean_1.fq.gz)" + "," + "\$(readlink -f ${id}_clean_2.fq.gz)"
    """

    merge_fastq_info.pl -s ${id} -a ${raw_info} -r ${rmhost_log} -o ${id}_reads_stat.txt

    echo ${id},${reads_path} > ${id}_clean_reads.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | sed 's/perl //g')
    END_VERSIONS
    
    """
    
}