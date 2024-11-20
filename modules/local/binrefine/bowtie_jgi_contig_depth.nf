
process BOWTIEJGICONTIGDEPTH {

    label 'process_single'

    input:
    path(index)
    path(reads)

    output:
    path("all_need_refine_bin_depth.xls"),emit:"binDepth"
    path("all_need_refine_bin_bin_bowtie2_log.txt"),emit:"bowtie2_log"
    path("all_need_refine_bin.sorted.bam"),emit:"sorted_bam"
    path("all_need_refine_bin_depth.list"),emit:"depth_list"

    when:
    task.ext.when == null || task.ext.when

    script:
    // def pre_command = params.single_end ? "zcat ${reads} > all_mash.fq" : "zcat ${reads[0]} ${reads[1]} > all_need_refine_bin.fq"
    // def reads_args = params.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def bowtie2_options = params.binset_profile_bowtie2_options ?: ""
    
    """
    
    zcat *bwa_mash_1.fq.gz > all_bwa_mash_1.fq

    zcat *bwa_mash_2.fq.gz > all_bwa_mash_2.fq

    bowtie2 -1 all_bwa_mash_1.fq -2 all_bwa_mash_2.fq -p ${task.cpus} -x index -S all_need_refine_bin.sam ${bowtie2_options} > all_need_refine_bin_bin_bowtie2_log.txt 2>&1
    
    samtools sort --threads ${task.cpus} all_need_refine_bin.sam --write-index -o all_need_refine_bin.sorted.bam
    
    jgi_summarize_bam_contig_depths --outputDepth all_need_refine_bin_depth.xls all_need_refine_bin.sorted.bam

    echo -e "all_need_refine_bin\t\$PWD/all_need_refine_bin_depth.xls" > all_need_refine_bin_depth.list

    #清理大文件
    #rm -rf *.sam *.bam
    #rm \$(readlink -f *.bt2)
    
    rm all_bwa_mash_*.fq

    """

}
