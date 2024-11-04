
include { BUILD4BIN } from '../../modules/local/binning/build4Bin'
include { BOWTIE2BIN } from '../../modules/local/binning/bowtie2Bin'
include { COVERM } from '../../modules/local/binning/coverm'
include { PIPELINEERROR } from '../../modules/local/common/pipeline_error'
include { PIPELINEEXIT } from '../../modules/local/common/pipeline_exit'


workflow BINABUNDANCE {
    
    take:
    sample_number
    clean_reads
    final_genomes       // RENAMEBIN.out.genomes
    bins_folder         // RENAMEBIN.out.folder
    method4coverm
    samplesheet
    
    main:

    if (!params.method4coverm) { exit 1, "The required parameter to calculate Bin Abundance is: --method4coverm" }

    //4. coverm calculates bin abundance.
    
    //bin taxonomic abundance estimation
 
    BUILD4BIN(final_genomes)

    BOWTIE2BIN(BUILD4BIN.out.bowtie2Index, clean_reads)
    bin_bowtie2 = BOWTIE2BIN.out.bowtie2_log.collect()
    depth_list = BOWTIE2BIN.out.depth_list
    bam_number = BOWTIE2BIN.out.sorted_bam.collect().flatten().filter { file -> !file.isEmpty() }.count()
    
    ch_coverm_method = Channel.of(method4coverm.split(","))

    COVERM(sample_number, bam_number, ch_coverm_method.combine(bins_folder), BOWTIE2BIN.out.sorted_bam.collect(), samplesheet)
    totalRelativeAbun = COVERM.out.abundance.collect().flatten().filter{ it.name.contains("relative_abundance") }
    totalCountAbun = COVERM.out.abundance.collect().flatten().filter{ it.name.contains("count") }
    totalMeanAbun = COVERM.out.abundance.collect().flatten().filter{ it.name.contains("trimmed_mean") }

    //Abnormal Execution.
    PIPELINEERROR("Bowtie2Bin", sample_number, bam_number)
    PIPELINEEXIT(PIPELINEERROR.out.log)

    binAundance_report = COVERM.out.report.collect().mix(PIPELINEERROR.out.log).collect()

    emit:
    bin_bowtie2
    totalRelativeAbun
    totalCountAbun
    totalMeanAbun
    binAundance_report
    depth_list

}