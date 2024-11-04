
include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'
include { BINFUNCTIONGENEID } from '../../modules/local/binning/bin_function_geneID'
include { BINSPECIFIEDFUNCTION as VFDBBINFUNCTION} from '../../modules/local/binning/bin_specified_function_geneID'
include { BINSPECIFIEDFUNCTION as CARDBINFUNCTION} from '../../modules/local/binning/bin_specified_function_geneID'
include { PIPELINEWARNING } from '../../modules/local/common/pipeline_warning'


workflow BINFUNCTION {
    take:
    gene_info         // channel: [ path(gene_info) ]
    clstr             // channel: [ path(clstr) ]
    annotations       // channel: [ path(All_annotations) ]
    final_genomes     // channel: RENAMEBIN.out.genomes
    vfdb_anno         // channel: GENESET.out.ch_vfdb_anno
    card_anno         // channel: GENESET.out.ch_card_anno

    main:

    /*
    * Verify the essential parameters for running this module
    */
    bin_function_essential_db = [params.cog_db_category, params.go_db_category, params.kegg_db_category, params.cazy_db_category]

    if(!checkEssentialParams(bin_function_essential_db)) { exit 1, "The required parameters to execute the Bin Function module are:\n --cog_db_category \n --go_db_category\n --kegg_db_category\n --cazy_db_category" }

    ch_cog_db_category = params2Channel(params.cog_db_category)
    ch_go_db_category = params2Channel(params.go_db_category)
    ch_kegg_db_category = params2Channel(params.kegg_db_category)
    ch_cazy_db_category =  params2Channel(params.cazy_db_category)


    // bin functional estimation 
    BINFUNCTIONGENEID (gene_info, clstr, annotations, final_genomes, ch_cog_db_category, ch_go_db_category, ch_kegg_db_category, ch_cazy_db_category)
    binFunction_geneID_report = BINFUNCTIONGENEID.out.eachBinFunction

    // binFunction with no result, output warning log
    PIPELINEWARNING("BinFunctionGeneID", BINFUNCTIONGENEID.out.warning)

    // input is optional
    VFDBBINFUNCTION(gene_info, clstr, vfdb_anno, final_genomes, "VFDB")
    CARDBINFUNCTION(gene_info, clstr, card_anno, final_genomes, "CARD")
    
    emit:
    binFunction_geneID_report

}