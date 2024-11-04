process BRADEPTH {
    
    tag "$binSample"
    
    label 'process_single'


    input:
    tuple val(binid),val(sampleid),val(binSample),path(reads1),path(reads2),path(bin_fa)


    output:
    tuple val(binid),val(sampleid),val(binSample),path("${binid}_${sampleid}.sorted.bam"),emit:"sorted_bam"
    tuple val(binid),val(sampleid),val(binSample),path("${binid}_${sampleid}_contig_depth.txt"),emit:"bin_depth"


    when:
    bin_fa.size() > 0

    script:
    def reads_args = params.single_end ? "-U ${reads1}" : "-1 ${reads1} -2 ${reads2}"
    def bowtie2_options = params.binning_bowtie2_options ?: ""
    
    """
    bowtie2-build ${bin_fa} ${binid} --threads ${task.cpus}
    bowtie2 ${reads_args} -p ${task.cpus} -x ${binid} -S ${binid}_${sampleid}.sam ${bowtie2_options} > ${binid}_${sampleid}_contig_bowtie2_log.txt 2>&1
    samtools sort --threads ${task.cpus} ${binid}_${sampleid}.sam --write-index -o ${binid}_${sampleid}.sorted.bam
    jgi_summarize_bam_contig_depths --outputDepth ${binid}_${sampleid}_contig_depth.txt ${binid}_${sampleid}.sorted.bam


    rm -rf *.sam *.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS

    """

}