process PREBINNING {
    
    tag "$id"
    
    label 'process_single'

    input:
    tuple val(id),path(contigs),path(reads)

    output:
    tuple val(id),path("${id}.sorted.bam"),emit:"sorted_bam"
    path("${id}.sorted.bam"),emit:"bam"
    tuple val(id),path("${id}.sorted.bam.csi"),emit:"sorted_bam_csi"
    tuple val(id),path("${id}_contig_depth.txt"),emit:"bin_depth"
    path("${id}_contig_bowtie2_log.txt"),emit:"bowtie2_log"


    when:
    contigs.size() > 0

    script:
    def reads_args = params.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def bowtie2_options = params.binning_bowtie2_options ?: ""
    
    """
    bowtie2-build ${contigs} ${id} --threads ${task.cpus}
    bowtie2 ${reads_args} -p ${task.cpus} -x ${id} -S ${id}.sam ${bowtie2_options} > ${id}_contig_bowtie2_log.txt 2>&1
    samtools sort --threads ${task.cpus} ${id}.sam --write-index -o ${id}.sorted.bam
    jgi_summarize_bam_contig_depths --outputDepth ${id}_contig_depth.txt ${id}.sorted.bam

    rm -rf *.sam *.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS

    """

    stub:
    """
    mkdir -p ${id}
    touch ${id}.sorted.bam
    touch ${id}.sorted.bam.csi
    touch ${id}_contig_depth.txt
    touch ${id}_contig_bowtie2_log.txt
    """

}