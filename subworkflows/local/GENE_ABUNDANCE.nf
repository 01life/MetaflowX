//
// gene abundance estimation
//

include { BUILD4GENE } from '../../modules/local/geneset/build4Gene'
include { BOWTIE2GENE } from '../../modules/local/geneset/bowtie2Gene'
include { MERGEGENEPROFILE } from '../../modules/local/geneset/merge_gene_profile'

include { PIPELINEERROR } from '../../modules/local/common/pipeline_error'
include { PIPELINEEXIT } from '../../modules/local/common/pipeline_exit'

workflow GENEABUNDANCE {
    take:
    sample_number
    clean_reads       // channel: [ val(id), [ reads1, reads2 ] ]
    samplesheet
    cdhit_cds
    gene_length

    main:

    BUILD4GENE(cdhit_cds)

    BOWTIE2GENE(BUILD4GENE.out.index, gene_length, clean_reads)
    geneset_bowtie2 = BOWTIE2GENE.out.bowtie2_log.collect()
    ch_gene_abun = BOWTIE2GENE.out.abundance.collect()
    finish_number_geneAbun = ch_gene_abun.flatten().filter { file -> !file.isEmpty() }.count()

    MERGEGENEPROFILE(sample_number, finish_number_geneAbun, ch_gene_abun, samplesheet)
    gene_profile = MERGEGENEPROFILE.out.merge
    abundance_report = MERGEGENEPROFILE.out.geneset_abundance_report

    // generate an error log and terminate the pipeline if GeneAbundanceCalculation error occurs
    PIPELINEERROR("GeneAbundanceCalculation", sample_number, finish_number_geneAbun)
    gene_abun_log = PIPELINEERROR.out.log
    PIPELINEEXIT(gene_abun_log)

    emit:
    gene_profile
    gene_abun_log
    abundance_report
    geneset_bowtie2

}

