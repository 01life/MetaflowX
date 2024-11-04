//
// geneset 
//

include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'
include { PRODIGAL } from '../../modules/local/geneset/prodigal'
// include { CDHIT } from '../../modules/local/geneset/cdhit'

include { GENEFILTER } from '../../modules/local/geneset/filter_gene_len'
include { SINGLECDHIT } from '../../modules/local/geneset/single_cdhit'
include { CDHITDIV } from '../../modules/local/geneset/multi_cdhit_1div'
include { ESTDIV } from '../../modules/local/geneset/multi_cdhit_2est_div'
include { CDHITEST2D } from '../../modules/local/geneset/multi_cdhit_3est2d'
include { CDHITCLSTR } from '../../modules/local/geneset/multi_cdhit_4create_mergeclstr'
include { MEGRCLSTR } from '../../modules/local/geneset/multi_cdhit_5mergeclstr'
include { PASTERESULT } from '../../modules/local/geneset/multi_cdhit_6paste_result'
include { CDHIT } from '../../modules/local/geneset/rename_cdhit'

include { EGGNOG } from '../../modules/local/geneset/eggnog'
include { BUILD4GENE } from '../../modules/local/geneset/build4Gene'
include { BOWTIE2GENE } from '../../modules/local/geneset/bowtie2Gene'
include { MERGEGENEPROFILE } from '../../modules/local/geneset/merge_gene_profile'
include { GENEPROFILE } from '../../modules/local/geneset/gene_profile'

include { RGI } from '../../modules/local/geneset/rgi'
include { ANTISMASH } from '../../modules/local/geneset/antismash'
include { CUSTOMNTDB } from '../../modules/local/geneset/custom_nucleotide_DB'
include { CUSTOMPRODB } from '../../modules/local/geneset/custom_protein_DB'
include { VFDB } from '../../modules/local/geneset/vfdb'
include { BIGMAP } from '../../modules/local/geneset/bigmap'
include { BIGMAPWARNING } from '../../modules/local/geneset/bigmap_warning'
include { MERGEBGCPROFILE } from '../../modules/local/geneset/merge_bgc_profile'
include { PIPELINEERROR as PIPELINEERROR_PRODIGAL } from '../../modules/local/common/pipeline_error'
include { PIPELINEEXIT as PIPELINEEXIT_PRODIGAL} from '../../modules/local/common/pipeline_exit'
include { PIPELINEERROR as PIPELINEERROR_GENEABUN} from '../../modules/local/common/pipeline_error'
include { PIPELINEEXIT as PIPELINEEXIT_GENEABUN} from '../../modules/local/common/pipeline_exit'
include { PIPELINEWARNING } from '../../modules/local/common/pipeline_warning'

workflow GENESET {
    take:
    sample_number
    contigs           // channel: [ val(id), path(contigs) ]
    clean_reads       // channel: [ val(id), [ reads1, reads2 ] ]
    samplesheet

    main:

    /*
    * Verify the essential parameters for running this module
    */
    geneset_essential_db = [params.eggnog_diamond_db, params.eggnog_mapper_db, params.cog_db_category, params.go_db_category, params.kegg_db_category, params.cazy_db_category]

    if(!checkEssentialParams(geneset_essential_db)) { exit 1, "The required parameters to execute the Geneset module are:\n --eggnog_diamond_db\n --eggnog_mapper_db\n --cog_db_category \n --go_db_category\n --kegg_db_category\n --cazy_db_category" }

    ch_eggnog_diamond_db = params2Channel(params.eggnog_diamond_db)
    ch_eggnog_mapper_db = params2Channel(params.eggnog_mapper_db)
    ch_cog_db_category = params2Channel(params.cog_db_category)
    ch_go_db_category = params2Channel(params.go_db_category)
    ch_kegg_db_category = params2Channel(params.kegg_db_category)
    ch_cazy_db_category =  params2Channel(params.cazy_db_category)


    ch_pipeline_log = Channel.empty()

    /* 
    * gene prediction
    */
    
    PRODIGAL(contigs)
    prodigal_faa = PRODIGAL.out.faa     // channel: [ val(id),path(faa) ]
    ch_prodigal_faa = PRODIGAL.out.noid_faa.collect()
    finish_number_prodigal = ch_prodigal_faa.flatten().filter { file -> !file.isEmpty() }.count()

    //CD-hit support multi tasks
    // CDHIT(sample_number, finish_number_prodigal, PRODIGAL.out.noid_cds.collect(), ch_prodigal_faa)
    GENEFILTER(sample_number, finish_number_prodigal, PRODIGAL.out.noid_cds.collect(), ch_prodigal_faa)

    // single cd-hit task
    SINGLECDHIT(GENEFILTER.out.allcds, GENEFILTER.out.allpep, GENEFILTER.out.single_task)

    // multi cd-hit task
    CDHITDIV(GENEFILTER.out.allcds, GENEFILTER.out.multi_task)
    ESTDIV(CDHITDIV.out.div.flatten())
    ch_div = CDHITDIV.out.div.flatten().map { file ->
            def fileName = file.name.replaceFirst(/^all.cds.fa.div-/, '')
                [fileName]
            }
    ch_est_input = ch_div.combine(CDHITDIV.out.div_tmp).combine(ESTDIV.out.div0.collect().toList())

    CDHITEST2D(ch_est_input)
    CDHITCLSTR(CDHITEST2D.out.clstr.collect(), CDHITEST2D.out.div_o.collect(), GENEFILTER.out.multi_task)
    ch_merge_input = CDHITCLSTR.out.clstr_order
        .splitText()
        .map { it -> 
            def in_list = it.trim().split("\t")[1]
            def clstr_num = it.trim().split("\t")[0]
            return [clstr_num,in_list]
        }
    MEGRCLSTR(ch_merge_input,CDHITEST2D.out.clstr.collect())
    PASTERESULT(MEGRCLSTR.out.subclstr.collect(), CDHITCLSTR.out.unique_fa, GENEFILTER.out.allpep)

    ch_cdhit_pep = SINGLECDHIT.out.pep.mix(PASTERESULT.out.pep).collect()
    ch_cdhit_cds = SINGLECDHIT.out.cds.mix(PASTERESULT.out.cds).collect()
    ch_cdhit_gene_length = SINGLECDHIT.out.gene_length.mix(PASTERESULT.out.gene_length).collect()
    ch_cdhit_gene_info = SINGLECDHIT.out.gene_info.mix(PASTERESULT.out.gene_info).collect()
    ch_cdhit_clstr = SINGLECDHIT.out.clstr.mix(PASTERESULT.out.clstr).collect()
    ch_cdhit_split = SINGLECDHIT.out.split.mix(PASTERESULT.out.split).collect()
    ch_cdhit_geneset_gene_report = SINGLECDHIT.out.geneset_gene_report.mix(PASTERESULT.out.geneset_gene_report).collect()

    CDHIT(ch_cdhit_pep, ch_cdhit_cds, ch_cdhit_gene_length, ch_cdhit_gene_info, ch_cdhit_clstr, ch_cdhit_split, ch_cdhit_geneset_gene_report)
    gene_info = CDHIT.out.gene_info     // channel: [ path(gene_info) ]
    clstr = CDHIT.out.clstr             // channel: [ path(clstr) ]
    gene_report = CDHIT.out.geneset_gene_report


    // generate an error log and terminate the pipeline if Prodigal error occurs
    PIPELINEERROR_PRODIGAL("Prodigal", sample_number, finish_number_prodigal)
    PIPELINEEXIT_PRODIGAL(PIPELINEERROR_PRODIGAL.out.log)
    
    /* 
    * gene functional annotation
    */
    eggNog_input = CDHIT.out.split.collect()
                            .flatten()
                            .map{ 
                                chunk -> [ chunk.baseName,chunk ]                                            
                            }
    
    EGGNOG(eggNog_input, ch_eggnog_diamond_db, ch_eggnog_mapper_db)
    
    /* 
    * gene abundance estimation
    */
    BUILD4GENE(CDHIT.out.cds)

    BOWTIE2GENE(BUILD4GENE.out.index, CDHIT.out.gene_length, clean_reads)
    geneset_bowtie2 = BOWTIE2GENE.out.bowtie2_log.collect()
    ch_gene_abun = BOWTIE2GENE.out.abundance.collect()
    finish_number_geneAbun = ch_gene_abun.flatten().filter { file -> !file.isEmpty() }.count()

    MERGEGENEPROFILE(sample_number, finish_number_geneAbun, ch_gene_abun, samplesheet)
    abundance_report = MERGEGENEPROFILE.out.geneset_abundance_report

    // generate an error log and terminate the pipeline if GeneAbundanceCalculation error occurs
    PIPELINEERROR_GENEABUN("GeneAbundanceCalculation", sample_number, finish_number_geneAbun)
    PIPELINEEXIT_GENEABUN(PIPELINEERROR_GENEABUN.out.log)

    ch_pipeline_log = PIPELINEERROR_PRODIGAL.out.log.mix(PIPELINEERROR_GENEABUN.out.log)

    /* 
    * gene functional abundance estimation
    */
    GENEPROFILE(EGGNOG.out.annotation.collect(), MERGEGENEPROFILE.out.merge, ch_cog_db_category, ch_go_db_category, ch_kegg_db_category, ch_cazy_db_category)
    geneset_total_abundance = GENEPROFILE.out.total_abundance
    annotation = GENEPROFILE.out.All_anontations // channel: [ path(All_annotations) ]
    function_report = GENEPROFILE.out.geneset_function_report
    ch_profile = GENEPROFILE.out.func_abun


    /* 
    * optional gene functional annotation
    */
    
    // VFDB annotation
    ch_vfdb_report = Channel.empty()
    vfdb_anno = Channel.empty()
    if(params.VFDB_db){
        ch_VFDB_db = params2Channel(params.VFDB_db)
        VFDB(CDHIT.out.pep, geneset_total_abundance, ch_VFDB_db)
        ch_vfdb_report = VFDB.out.vfdb_info
        vfdb_anno = VFDB.out.vfdb_anno
    }
    
    // RGI annotation
    ch_card_report = Channel.empty()
    card_anno = Channel.empty()
    if(params.CARD_db){
        ch_CARD_db = params2Channel(params.CARD_db)
        RGI(CDHIT.out.pep, geneset_total_abundance, ch_CARD_db)
        ch_card_report = RGI.out.rgi_info
        card_anno = RGI.out.rgi_anno
    }

    // antiSMASH annotation
    ch_bigmap_report = Channel.empty()
    if(params.bigspace_db){
        ch_bigspace_db = params2Channel(params.bigspace_db)

        ANTISMASH(contigs)
        
        BIGMAP(clean_reads.join(ANTISMASH.out.all_bin), ch_bigspace_db)
        ch_bigmap_out = BIGMAP.out.corecov.collect()
        finish_number_bigmap = ch_bigmap_out.flatten().filter { file -> !file.isEmpty() }.count()

        // generate a warning log in the event of a BIGMAP error
        BIGMAPWARNING(sample_number, finish_number_bigmap)
        PIPELINEWARNING("BiG-MAP", BIGMAPWARNING.out.log)
        ch_pipeline_log = ch_pipeline_log.mix(BIGMAPWARNING.out.log)

        MERGEBGCPROFILE(BIGMAP.out.BGCfamily.collect(), samplesheet)
        ch_bigmap_report = MERGEBGCPROFILE.out.bigmap_info

    }
    
    // cutsom nucleotide database annotation
    ch_customNT_report = Channel.empty()
    if(params.nucleotide_db){
        ch_nucleotide_db = params2Channel(params.nucleotide_db)
        CUSTOMNTDB(CDHIT.out.cds, geneset_total_abundance, ch_nucleotide_db)
        ch_customNT_report = CUSTOMNTDB.out.customnt_info
    }

    // custom protein database annotation
    ch_customPR_report = Channel.empty()
    if(params.protein_db){
        ch_protein_db = params2Channel(params.protein_db)
        CUSTOMPRODB(CDHIT.out.pep, geneset_total_abundance, ch_protein_db)
        ch_customPR_report = CUSTOMPRODB.out.custompr_info
    }
    
    optional_report = ch_vfdb_report.mix(ch_card_report, ch_bigmap_report, ch_customNT_report, ch_customPR_report)
    
    geneset_report = gene_report.mix(abundance_report, function_report, optional_report, ch_pipeline_log).collect()

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

