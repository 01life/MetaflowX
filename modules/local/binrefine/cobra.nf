process COBRA {

    tag "$bin_id"
    
    label 'process_single'

    publishDir "${params.outdir}/08.BinOptimization/081.BinRefine/COBRA/",mode:'copy'

    input: 
    tuple val(bin_id),path(bin_fa),path(mash_reads1),path(mash_reads2),path(allcontigs),path(index)

    output:
    path("COBRA_${bin_id}.fa"),emit:"bin_cobra"

    script:
    def bowtie2_options = params.binning_bowtie2_options ?: ""
    def assembler = params.metaspades ? "metaspades" : "megahit"
    def cobra_options = params.cobra_options ?: ""
    
    """
    # step1 obtain the coverage file

    bowtie2 -1 ${mash_reads1} -2 ${mash_reads2} -p ${task.cpus} -x index -S ${bin_id}.sam ${bowtie2_options} > ${bin_id}_contig_bowtie2_log.txt 2>&1

    samtools sort --threads ${task.cpus} ${bin_id}.sam --write-index -o ${bin_id}.sorted.bam

    jgi_summarize_bam_contig_depths --outputDepth ${bin_id}_original_contig_depth.txt ${bin_id}.sorted.bam

    coverage.transfer.py -i ${bin_id}_original_contig_depth.txt -o ${bin_id}_contig_depth.txt


    # step2 run cobra

    cobra-meta -f ${allcontigs} -q ${bin_fa} -o ${bin_id}.COBRA -c ${bin_id}_contig_depth.txt -m ${bin_id}.sam -a ${assembler} -t ${task.cpus} ${cobra_options}

    # Check if contig.new.fa exists.
    if [ ! -f "${bin_id}.COBRA/contig.new.fa" ]; then
        touch "COBRA_${bin_id}.fa"
    else
        cp ${bin_id}.COBRA/contig.new.fa COBRA_${bin_id}.fa
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra //; s/ .*\$//')
    END_VERSIONS

    """

}