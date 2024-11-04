
process BOWTIE2BIN {
    
    tag "$id"

    label 'process_single'

    input:
    path(index)
    tuple val(id),path(reads)

    output:
    tuple val(id),path("${id}_depth.xls"),emit:"binDepth"
    path("${id}_bin_bowtie2_log.txt"),emit:"bowtie2_log"
    path("${id}.sorted.bam"),emit:"sorted_bam"
    path("${id}_depth.list"),emit:"depth_list"

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_args = params.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def bowtie2_options = params.binset_profile_bowtie2_options ?: ""
    
    """

    bowtie2 ${reads_args} -p ${task.cpus} -x all.bin.fa -S ${id}.sam ${bowtie2_options} > ${id}_bin_bowtie2_log.txt 2>&1
    
    samtools sort --threads ${task.cpus} ${id}.sam --write-index -o ${id}.sorted.bam
    
    jgi_summarize_bam_contig_depths --outputDepth ${id}_depth.xls ${id}.sorted.bam

    echo -e "${id}\t\$PWD/${id}_depth.xls" > ${id}_depth.list

    #Clean up intermediate files.
    #rm -rf *.sam *.bam
    #rm \$(readlink -f *.bt2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS

    """

}
