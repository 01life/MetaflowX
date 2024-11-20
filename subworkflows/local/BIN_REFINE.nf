include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'

include { GETBINMAPREADS } from '../../modules/local/binrefine/get_bin_map_reads'
include { DEEPURIFYCLEAN } from '../../modules/local/binrefine/deepurify_clean'
include { POSTDEEPURIFYCLEAN } from '../../modules/local/binrefine/post_deepurify_clean'
include { DEEPURIFYREBIN } from '../../modules/local/binrefine/deepurify_rebin'
include { COBRA } from '../../modules/local/binrefine/cobra'
include { GETBINMASHFQ } from '../../modules/local/binrefine/get_bin_mash_fq'
include { REBINSTAT } from '../../modules/local/binrefine/rebin_stat'
include { BUILD4GENE as BUILD4ORGBINS } from '../../modules/local/binrefine/build4Gene'
include { BOWTIEJGICONTIGDEPTH } from '../../modules/local/binrefine/bowtie_jgi_contig_depth'
include { BUILD4GENE } from '../../modules/local/binrefine/build4Gene'
include { MERGEREFINEBIN } from '../../modules/local/binrefine/merge_refine_bin'
include { CHECKM2 } from '../../modules/local/binning/checkm2'
include { FILTERBINS } from '../../modules/local/binning/filter_bins'
include { PIPELINEWARNING as CHECKM2_WARNING } from '../../modules/local/common/pipeline_warning'
include { FILTERREFINE } from '../../modules/local/binrefine/refine_bin_filter'


workflow BINREFINE {
    take:
    clean_reads
    contigs //ch_contigs
    filter_bin
    filter_bin_info
    filter_bin_mash_fq

    main:
    
    /*
    * Verify the essential parameters for running this module
    */
    binner_essential_db = [params.checkm2_db]

    if(!checkEssentialParams(binner_essential_db)) { exit 1, "The required parameter to execute the Bin Refine module is: checkm2_db" }

    ch_checkm2_db = params2Channel(params.checkm2_db)


    // GETBINMAPREADS(xxx, xxx) // [bin_id, bin_fa, ref_fa, ha_fq_id, fq_id_list], ch_clean_reads.collect()

    ch_data = filter_bin_mash_fq.splitText(keepHeader: true)
        .map { it -> 
            def split = it.split("\t")
            def bin_id = split[0].trim()
            def bin_fa = split[1].trim()
            def ref_fa = split[2].trim()
            def ha_fq_id = split[3].trim()
            def fq_id_list = split[4].trim()
            return [ bin_id, bin_fa, ref_fa, ha_fq_id, fq_id_list ]
        }

    ch_clean_fq = clean_reads.map{ it -> it[1] }.collect()

    GETBINMAPREADS(ch_clean_fq,ch_data)


    ch_deepurify_res = Channel.empty()
    
    if(params.deepurify_module == "re-bin"){
       
        all_org_bins = ch_data.map{ it -> return file(it[1])}.collectFile(name:"all_org_bins.fa")

        BUILD4ORGBINS(all_org_bins)

        BOWTIEJGICONTIGDEPTH(BUILD4ORGBINS.out.index, GETBINMAPREADS.out.map.map{ it -> [it[1], it[2]] }.flatten().collect())

        DEEPURIFYREBIN(all_org_bins, BOWTIEJGICONTIGDEPTH.out.sorted_bam, ch_checkm2_db)
        ch_deepurify_res = DEEPURIFYREBIN.out.res

    }else if(params.deepurify_module == "clean"){
        
        DEEPURIFYCLEAN(filter_bin, ch_checkm2_db)
        ch_deepurify_res = DEEPURIFYCLEAN.out.res

    }else{
        exit 1, "deepurify: error: argument command: invalid choice: 'rebin' (choose from 'clean', 're-bin') !"
    }

    // handle depurify result
    POSTDEEPURIFYCLEAN(filter_bin, filter_bin_mash_fq, ch_deepurify_res)

    ch_deepurify_bin = POSTDEEPURIFYCLEAN.out.deepurify_bin
    ch_deepurify_bin_mash_fq = POSTDEEPURIFYCLEAN.out.deepurify_bin_mash_fq


    // Merge contig information of all samples.
    all_contigs = contigs.map{ _,file -> file}.collectFile(name:"all_contigs.fa")

    BUILD4GENE(all_contigs)

    ch_cobra_in = ch_deepurify_bin.flatten()
        .map{ it -> 
            def fname = file(it).name
            def id = fname.replace("Deepurify_", "").replace(".fa", "")
            return [id, it]
        }
        .join(GETBINMAPREADS.out.map)
        .combine(all_contigs)
        .combine(BUILD4GENE.out.index.toList())

    COBRA(ch_cobra_in) 

    MERGEREFINEBIN(ch_deepurify_bin, COBRA.out.bin_cobra.collect())

    after_refine_bin_mash_fq = Channel.empty()

    ch_checkm2_in = MERGEREFINEBIN.out.refine_bin
    
    //Quality Control.
    CHECKM2(ch_checkm2_in, ch_checkm2_db)
    CHECKM2_WARNING("CheckM2", CHECKM2.out.checkm2_error.collect())

    //Filtering.
    filter_bins_input = ch_checkm2_in.join(CHECKM2.out.quality_report,by:0)
    FILTERBINS (filter_bins_input)
    BIN_FILTER_WARNING("COBRA_checkm2", FILTERBINS.out.filtered_log.collect())

 
    FILTERREFINE(filter_bin_info, FILTERBINS.out.QS_quality_report, MERGEREFINEBIN.out.refine_bin_list, ch_deepurify_bin_mash_fq) 
    
    after_refine_bin_mash_fq = FILTERREFINE.out.after_refine_bin_mash_fq


    emit:
    after_refine_bin_mash_fq
}