//
// geneset analysis
//

include { GENEPREDICTION } from './GENE_PREDICTION'
include { GENESETCONSTRUCTION } from './GENESET_CONSTRUCTION'
include { GENEABUNDANCE } from './GENE_ABUNDANCE'
include { GENEFUNCTION } from './GENE_FUNCTION'


workflow GENESET {
    take:
    sample_number
    contigs           // channel: [ val(id), path(contigs) ]
    clean_reads       // channel: [ val(id), [ reads1, reads2 ] ]
    samplesheet     

    main:

    geneset_report = Channel.empty()
    geneset_log = Channel.empty()

    // default execute gene prediction
    GENEPREDICTION(sample_number, contigs)
    prodigal_faa = GENEPREDICTION.out.prodigal_faa
    prodigal_cds = GENEPREDICTION.out.prodigal_cds
    // prodigal_noid_faa = GENEPREDICTION.out.prodigal_noid_faa
    // prodigal_noid_cds = GENEPREDICTION.out.prodigal_noid_cds
    prodigal_log = GENEPREDICTION.out.prodigal_log

    clstr = Channel.empty()
    split_fa = Channel.empty()
    cds = Channel.empty()
    pep = Channel.empty()
    gene_info = Channel.empty()
    gene_length = Channel.empty()
    if(params.geneset_construction){
        GENESETCONSTRUCTION(sample_number, prodigal_cds, prodigal_faa)
        clstr = GENESETCONSTRUCTION.out.clstr
        split_fa = GENESETCONSTRUCTION.out.split_fa
        cds = GENESETCONSTRUCTION.out.cds
        pep = GENESETCONSTRUCTION.out.pep
        gene_info = GENESETCONSTRUCTION.out.gene_info
        gene_length = GENESETCONSTRUCTION.out.gene_length
        gene_report = GENESETCONSTRUCTION.out.gene_report
        geneset_report = geneset_report.mix(gene_report)
    }

    gene_profile = Channel.empty()
    geneset_bowtie2 = Channel.empty()
    abundance_report = Channel.empty()
    gene_abun_log = Channel.empty()
    if(params.gene_abundance_calculation){
        GENEABUNDANCE(sample_number, clean_reads, samplesheet, cds, gene_length)
        gene_profile = GENEABUNDANCE.out.gene_profile
        geneset_bowtie2 = GENEABUNDANCE.out.geneset_bowtie2
        abundance_report = GENEABUNDANCE.out.abundance_report
        gene_abun_log = GENEABUNDANCE.out.gene_abun_log
    }

    ch_profile = Channel.empty()
    function_report = Channel.empty()
    annotation = Channel.empty()
    vfdb_anno = Channel.empty()
    card_anno = Channel.empty()
    optional_report = Channel.empty()
    bigmap_warning = Channel.empty()
    if(params.gene_function){
        GENEFUNCTION(sample_number, contigs, clean_reads, samplesheet, cds, pep, split_fa, gene_profile)
        ch_profile = GENEFUNCTION.out.function_abundance
        function_report = GENEFUNCTION.out.function_report
        annotation = GENEFUNCTION.out.eggnog_annotation
        vfdb_anno = GENEFUNCTION.out.vfdb_anno
        card_anno = GENEFUNCTION.out.card_anno
        optional_report = GENEFUNCTION.out.optional_report
        bigmap_warning = GENEFUNCTION.out.bigmap_warning
    }
    
    geneset_log = prodigal_log.mix(gene_abun_log, bigmap_warning)
    geneset_report = geneset_report.mix(abundance_report, function_report, optional_report, geneset_log).collect()

    emit:
    gene_info 
    clstr
    prodigal_faa
    annotation
    geneset_bowtie2
    geneset_report
    ch_profile
    vfdb_anno //VFDB_annotation.xls
    card_anno //CARD_annotation.xls
    
}

