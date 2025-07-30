//
// geneset functional analysis
//

include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'

include { EGGNOG } from '../../modules/local/geneset/eggnog'
include { GENEPROFILE } from '../../modules/local/geneset/gene_profile'

include { RGI } from '../../modules/local/geneset/rgi'
include { ANTISMASH } from '../../modules/local/geneset/antismash'
include { CUSTOMNTDB } from '../../modules/local/geneset/custom_nucleotide_DB'
include { CUSTOMPRODB } from '../../modules/local/geneset/custom_protein_DB'
include { VFDB } from '../../modules/local/geneset/vfdb'
include { BIGMAP } from '../../modules/local/geneset/bigmap'
include { BIGMAPWARNING } from '../../modules/local/geneset/bigmap_warning'
include { MERGEBGCPROFILE } from '../../modules/local/geneset/merge_bgc_profile'
include { PIPELINEWARNING } from '../../modules/local/common/pipeline_warning'

workflow GENEFUNCTION {
    take:
    sample_number
    contigs           // channel: [ val(id), path(contigs) ]
    clean_reads       // channel: [ val(id), [ reads1, reads2 ] ]
    samplesheet
    cdhit_cds
    cdhit_pep         // CDHIT.out.pep
    cdhit_split_fa    // CDHIT.out.split.collect() 
    gene_profile      // MERGEGENEPROFILE.out.merge


    main:
    /*
    * Verify the essential parameters for running this module
    */
    gene_function_essential_db = [params.eggnog_diamond_db, params.eggnog_mapper_db, params.cog_db_category, params.go_db_category, params.kegg_db_category, params.cazy_db_category]

    if(!checkEssentialParams(gene_function_essential_db)) { exit 1, "The required parameters to execute the GENEFUNCTION module are:\n --eggnog_diamond_db\n --eggnog_mapper_db\n --cog_db_category \n --go_db_category\n --kegg_db_category\n --cazy_db_category" }

    ch_eggnog_diamond_db = params2Channel(params.eggnog_diamond_db)
    ch_eggnog_mapper_db = params2Channel(params.eggnog_mapper_db)
    ch_cog_db_category = params2Channel(params.cog_db_category)
    ch_go_db_category = params2Channel(params.go_db_category)
    ch_kegg_db_category = params2Channel(params.kegg_db_category)
    ch_cazy_db_category =  params2Channel(params.cazy_db_category)

    /* 
    * gene functional annotation
    */
    eggNog_input = cdhit_split_fa.flatten()
                    .map{ 
                        chunk -> [ chunk.baseName, chunk]                                            
                    }

    EGGNOG(eggNog_input, ch_eggnog_diamond_db, ch_eggnog_mapper_db)
    
    /* 
    * get special function database profile from eggNOG-mapper result
    */
    GENEPROFILE(EGGNOG.out.annotation.collect(), gene_profile, ch_cog_db_category, ch_go_db_category, ch_kegg_db_category, ch_cazy_db_category)
    geneset_total_abundance = GENEPROFILE.out.total_abundance
    eggnog_annotation = GENEPROFILE.out.All_anontations // channel: [ path(All_annotations) ]
    function_abundance = GENEPROFILE.out.func_abun
    function_report = GENEPROFILE.out.geneset_function_report


    /* 
    * optional gene functional annotation
    */
    
    // VFDB annotation
    ch_vfdb_report = Channel.empty()
    vfdb_anno = Channel.empty()
    if(params.VFDB_db){
        ch_VFDB_db = params2Channel(params.VFDB_db)
        VFDB(cdhit_pep, geneset_total_abundance, ch_VFDB_db)
        ch_vfdb_report = VFDB.out.vfdb_info
        vfdb_anno = VFDB.out.vfdb_anno
    }
    
    // RGI annotation
    ch_card_report = Channel.empty()
    card_anno = Channel.empty()
    if(params.CARD_db){
        ch_CARD_db = params2Channel(params.CARD_db)
        RGI(cdhit_pep, geneset_total_abundance, ch_CARD_db)
        ch_card_report = RGI.out.rgi_info
        card_anno = RGI.out.rgi_anno
    }

    // antiSMASH annotation
    ch_bigmap_report = Channel.empty()
    bigmap_warning = Channel.empty()
    if(params.bigspace_db){
        ch_bigspace_db = params2Channel(params.bigspace_db)

        ANTISMASH(contigs)
        
        BIGMAP(clean_reads.join(ANTISMASH.out.all_bin), ch_bigspace_db)
        ch_bigmap_out = BIGMAP.out.corecov.collect()
        finish_number_bigmap = ch_bigmap_out.flatten().filter { file -> !file.isEmpty() }.count()

        // generate a warning log in the event of a BIGMAP error
        BIGMAPWARNING(sample_number, finish_number_bigmap)
        PIPELINEWARNING("BiG-MAP", BIGMAPWARNING.out.log)
        bigmap_warning = BIGMAPWARNING.out.log

        MERGEBGCPROFILE(BIGMAP.out.BGCfamily.collect(), samplesheet)
        ch_bigmap_report = MERGEBGCPROFILE.out.bigmap_info

    }
    
    // cutsom nucleotide database annotation
    ch_customNT_report = Channel.empty()
    if(params.nucleotide_db){
        ch_nucleotide_db = params2Channel(params.nucleotide_db)
        CUSTOMNTDB(cdhit_cds, geneset_total_abundance, ch_nucleotide_db)
        ch_customNT_report = CUSTOMNTDB.out.customnt_info
    }

    // custom protein database annotation
    ch_customPR_report = Channel.empty()
    if(params.protein_db){
        ch_protein_db = params2Channel(params.protein_db)
        CUSTOMPRODB(cdhit_pep, geneset_total_abundance, ch_protein_db)
        ch_customPR_report = CUSTOMPRODB.out.custompr_info
    }

    optional_report = ch_vfdb_report.mix(ch_card_report, ch_bigmap_report, ch_customNT_report, ch_customPR_report)

    emit:
    function_abundance
    function_report
    eggnog_annotation
    vfdb_anno //VFDB_annotation.xls
    card_anno //CARD_annotation.xls
    bigmap_warning
    optional_report
    
}

