//
// BUILD BIN SET
//
include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'
include { MAXBIN2 } from '../../modules/local/binning/maxbin2'
include { METABAT2 } from '../../modules/local/binning/metabat2'
include { CONCOCT } from '../../modules/local/binning/concoct'
include { SEMIBIN2 } from '../../modules/local/binning/semibin2'
include { METABINNER } from '../../modules/local/binning/metabinner'
include { BINNY } from '../../modules/local/binning/binny'
include { COMEBIN } from '../../modules/local/binning/comebin'
include { VAMBBIN } from '../../modules/local/binning/vamb_bin'

include { CONTIGS_TAXONOMY } from './CONTIGS_TAXONOMY'

include { MULTIBINNERWARNING } from '../../modules/local/binning/multi_binner_warning'
include { COMBINEBINNER } from '../../modules/local/binning/combine_binner_contig2tsv'
include { RENAMECHECKM2 } from '../../modules/local/binning/renamecheckm2'

include { CHECKM2 as COMBINECHECKM2 } from '../../modules/local/binning/checkm2'
include { SELECTPERMUTATION } from '../../modules/local/binning/select_bin_permutation'

include { DASTOOL } from '../../modules/local/binning/das_tool'
include { MAGSCOT } from '../../modules/local/binning/magscot'


include { GETSAMPLEBINMAP } from '../../modules/local/binning/get_sample_bin_map'
include { CHECKM2 } from '../../modules/local/binning/checkm2'
include { QUAST_SAMPLE_BINS } from '../../modules/local/binning/quast_sample_bins'
include { PICKPBODBO } from '../../modules/local/binning/pick_best_bin'
include { FILTERBINS } from '../../modules/local/binning/filter_bins'
include { PREBINNING } from '../../modules/local/binning/prebinning'
include { POSTBINNING } from '../../modules/local/binning/postbinning'
include { DREP } from '../../modules/local/binning/drep'
include { GALAH } from '../../modules/local/binning/galah'
include { RENAMEBIN } from '../../modules/local/binning/rename_bin'
include { PIPELINEWARNING as BINNER_WARNING } from '../../modules/local/common/pipeline_warning'
include { PIPELINEWARNING as DASTOOL_WARNING } from '../../modules/local/common/pipeline_warning'
include { PIPELINEWARNING as CHECKM2_WARNING } from '../../modules/local/common/pipeline_warning'
include { PIPELINEWARNING as BIN_FILTER_WARNING } from '../../modules/local/common/pipeline_warning'
include { NOBINSWARNING } from '../../modules/local/binning/no_bins_warning'


workflow BINNER {
    take:
    contigs           // channel: [ val(id), [ contigs ] ]
    clean_reads       // channel: [ val(id), [ reads1, read2 ] ]
    prodigal_faa      // channel: [ val(id), path(faa) ]
    contig_taxonomy   // channel: [ val(id), path(taxonomy) ]
    contig_map        // channel: [ val(id), path(contig_map) ]

    main:

    /*
    * Verify the essential parameters for running this module
    */
    binner_essential_db = [params.checkm2_db]

    if(!checkEssentialParams(binner_essential_db)) { exit 1, "The required parameter to execute the Binner module is: checkm2_db" }

    ch_checkm2_db = params2Channel(params.checkm2_db)
    ch_magscot_folder = params2Channel(params.magscot_folder)

    contig_map.view()


    /*
    * binning 
    */
    
    def ch_binner_list = [params.maxbin2, params.metabat2, params.concoct, params.metabinner, params.semibin2, params.binny, params.comebin, params.vamb].findAll{ it==true }

    def binner_number = ch_binner_list.size()

    println(binner_number)

    def single_binner_result = Channel.empty()

    def multi_binner_log = Channel.empty()

    //Exit the program if no binning tool is specified.
    if ( binner_number == 0 ){
        exit 1, 'Error: tools for binning is not set !'
    }

    bin_input = contigs.join(clean_reads)

    concoct_input = contigs
    metabat_input = contigs

    ID = contigs.map{it -> it[0]}.first() 

    def preFlag = [params.metabat2, params.concoct, params.metabinner, params.semibin2, params.binny, params.comebin, params.vamb].findAll{ it==true }

    //For binning tools other than maxbin, first execute prebinning.
    if( preFlag.size() > 0 ){
        PREBINNING(bin_input)
        contig_bowtie2 = PREBINNING.out.bowtie2_log.collect()
        concoct_input =  contigs.join(PREBINNING.out.sorted_bam).join(PREBINNING.out.sorted_bam_csi)
        metabat_input =  contigs.join(PREBINNING.out.bin_depth)
        metabinner_input = metabat_input
        vamb_input = metabat_input
        semibin2_input = contigs.join(PREBINNING.out.sorted_bam)
        binny_input = semibin2_input
        comebin_input = semibin2_input
    }

    if ( params.maxbin2 ){
        MAXBIN2 (bin_input)
        single_binner_result = MAXBIN2.out.bins
    }

    if ( params.metabat2 ){
        METABAT2 (metabat_input)
        single_binner_result = METABAT2.out.bins
    }

    if ( params.concoct ){
        CONCOCT (concoct_input)
        single_binner_result = CONCOCT.out.bins
    }
    
    //Add new binning tools: metabinner, semibin2, binny, comebin.
    if ( params.metabinner ){
        METABINNER (metabinner_input)
        single_binner_result = METABINNER.out.bins
    }

    if ( params.semibin2 ){
        SEMIBIN2 (semibin2_input)
        single_binner_result = SEMIBIN2.out.bins
    }

    if ( params.binny ){
        BINNY (binny_input)
        single_binner_result = BINNY.out.bins
    }

    if( params.comebin){
        COMEBIN(comebin_input)
        single_binner_result = COMEBIN.out.bins
    }

    if( params.vamb){
        VAMBBIN(vamb_input)
        single_binner_result = VAMBBIN.out.bins
    }
    
    def ch_maxbin2_tsv = params.maxbin2 ?  MAXBIN2.out.tsv : channel.empty()
    def ch_metabat2_tsv = params.metabat2 ? METABAT2.out.tsv : channel.empty()
    def ch_concoct_tsv = params.concoct ? CONCOCT.out.tsv : channel.empty()
    def ch_metabinner_tsv = params.metabinner ?  METABINNER.out.tsv : channel.empty()
    def ch_semibin2_tsv = params.semibin2 ? SEMIBIN2.out.tsv : channel.empty()
    def ch_binny_tsv = params.binny ? BINNY.out.tsv : channel.empty()
    def ch_comebin_tsv = params.comebin ? COMEBIN.out.tsv : channel.empty()
    def ch_vamb_tsv = params.vamb ? VAMBBIN.out.tsv : channel.empty()

    // BinsContigs
    def ch_maxbin2_bc = params.maxbin2 ?  MAXBIN2.out.BinsContigs : channel.empty()
    def ch_metabat2_bc = params.metabat2 ? METABAT2.out.BinsContigs : channel.empty()
    def ch_concoct_bc = params.concoct ? CONCOCT.out.BinsContigs : channel.empty()
    def ch_metabinner_bc = params.metabinner ?  METABINNER.out.BinsContigs : channel.empty()
    def ch_semibin2_bc = params.semibin2 ? SEMIBIN2.out.BinsContigs : channel.empty()
    def ch_binny_bc = params.binny ? BINNY.out.BinsContigs : channel.empty()
    def ch_comebin_bc = params.comebin ? COMEBIN.out.BinsContigs : channel.empty()
    def ch_vamb_bc = params.vamb ? VAMBBIN.out.BinsContigs : channel.empty()



    //*****************
    //Select multiple binning tools.
    //*****************

    //  two way to optimize Bins
    ch_best_bin     = Channel.empty()
    ch_DBO_qs       = Channel.empty()
    ch_PBO_qs       = Channel.empty()
    ch_DBO_bin      = Channel.empty()
    ch_PBO_bin      = Channel.empty()
    ch_dastool_in   = Channel.empty()
    ch_magscot_in   = Channel.empty()
    combine_bin_ts  = Channel.empty()
    ch_checkm2_in   = Channel.empty()

    def optimizeBins_method = ''

    if (binner_number>1){
        //Merge TSV: [id, [tsv]].
        def bins_tsv = ch_maxbin2_tsv.mix(ch_metabat2_tsv, ch_concoct_tsv, ch_metabinner_tsv, ch_semibin2_tsv, ch_binny_tsv, ch_comebin_tsv, ch_vamb_tsv).groupTuple(by: 0)

        //Merge BinsContigs: [id, [tsv]].
        def bins_bc = ch_maxbin2_bc.mix(ch_metabat2_bc, ch_concoct_bc, ch_metabinner_bc, ch_semibin2_bc, ch_binny_bc, ch_comebin_bc, ch_vamb_bc).groupTuple(by: 0)

        //When multiple binning tools encounter an exception, output a log reminder.
        binner_tsv_number = bins_tsv.map { id, files ->   
            // Find non-empty files (i.e., binning is normal).
            def nonEmptyFiles = files.findAll { it ->  !it.isEmpty() }  
            // Return sample ID, number of binning tools, and number of binning result TSV files.
            [id, binner_number, nonEmptyFiles.size()]  
        }
        // Execute only when the number of binning tools != number of binning results.
        MULTIBINNERWARNING(binner_tsv_number)
        BINNER_WARNING("Binning", MULTIBINNERWARNING.out.log.collect())

        // *two way to optimize the best bins* //
        // DASTool-based Bin Optimizer
        if ( params.DASToolBinOptimizer  ){
            optimizeBins_method = 'DBO'
            //dastools input [id, contigs, faa, [binner tsv]]
            ch_dastool_in = contigs.join(prodigal_faa).join(bins_tsv)
            ch_magscot_in = contigs.join(prodigal_faa).join(bins_bc)
        }

        // Permutation-based Bin Optimizer
        ch_PBO_tsv = Channel.empty()
        if ( params.PermutationBinOptimizer ){
            optimizeBins_method = 'PBO'

            // combine binner tsv [id, [tsv]]
            ch_combinebin_in = bins_tsv.join(PREBINNING.out.bin_depth).join(contigs)
            COMBINEBINNER(ch_combinebin_in)
            combine_bin_tsv = COMBINEBINNER.out.combine_bin_tsv

            // checkm2 input , create unique id [id,floder]
            // ch_combinecheck_in = COMBINEBINNER.out.combine_bin_fa
            //                         .flatMap { id, folders ->
            //                             folders.collect { folder ->
            //                                 // Extract the folder name from the folder path
            //                                 def foldername = folder.getName()
            //                                 // Concatenate id with foldername and return the new tuple
            //                                 tuple("${id}-${foldername}", folder)
            //                             }
            //                         }

            ch_combinecheck_in = COMBINEBINNER.out.combine_bin_fa
                .flatMap { id, folders ->
                    // Handle multiple folders or single folder
                    if (folders instanceof List) {
                        // If folders is a list, process each folder
                        return folders.collect { folder ->
                            // Extract name from folder path
                            def foldername = folder.getName()
                            // Combine id and folder name, return new tuple
                            tuple("${id}-${foldername}", folder)
                        }
                    } else {
                        // If there's only one folder, process it directly
                        def foldername = folders.getName()
                        return [tuple("${id}-${foldername}", folders)]
                    }
                }

            
            COMBINECHECKM2(ch_combinecheck_in, ch_checkm2_db)

            RENAMECHECKM2(COMBINECHECKM2.out.quality_report)

            //[id, [qs list]]
            ch_combin_check_report = RENAMECHECKM2.out.new_quality_report
                .map { id, filepath ->
                    // Extract the base identifier by splitting at the first underscore
                    def base_id = id.split('-')[0]
                    // Return the new tuple with the base_id and the original filepath
                    tuple(base_id, filepath)
                    }.groupTuple(by: 0)
            ch_permution_in = COMBINEBINNER.out.allcontigs2bin.join( ch_combin_check_report ).join(contigs)

            // pcik the best combination [id,[tsvs],[qs reports],contig]
            SELECTPERMUTATION(ch_permution_in)
            ch_PBO_tsv = SELECTPERMUTATION.out.contigs2bin
            ch_PBO_qs = SELECTPERMUTATION.out.bin_qs
            ch_PBO_bin = SELECTPERMUTATION.out.best_bin
            ch_PBO_bc = SELECTPERMUTATION.out.BinsContigs



        }

        if ( params.ContigTaxonomyOptimizer ){

            CONTIGS_TAXONOMY(contigs,PREBINNING.out.bin_depth,clean_reads,contig_taxonomy,contig_map)

        }


        if (  params.DASToolBinOptimizer && params.PermutationBinOptimizer ){
            optimizeBins_method = 'DBO-PBO'
            //newtsv = [[binner tsvs] ,  [combine tsvs] ]
            //dastools input [id, contigs, faa, [newtsvs]]
            // new_tsv = bins_tsv.mix(ch_PBO_tsv) //.groupTuple(by: 0)
            new_tsv = bins_tsv.join(ch_PBO_tsv)
                .map { item ->  
                    def id = item[0]  
                    def tsv = item[1]  
                    def tsv2 = item[2]  
                    // Add values to the list.
                    tsv.add(tsv2)  
                    // Return new format.
                    return [id, tsv]  
                } 
            ch_dastool_in = contigs.join(prodigal_faa).join(new_tsv)

            new_bc = bins_bc.join(ch_PBO_bc)
                .map { item ->  
                    def id = item[0]  
                    def bc = item[1]  
                    def bc2 = item[2]  
                    // Add values to the list.
                    bc.add(bc2)  
                    // Return new format.
                    return [id, bc]  
                }
            ch_magscot_in = contigs.join(prodigal_faa).join(new_bc)

        }
        DASTOOL(ch_dastool_in)
        //Output log for dastool exceptions.
        DASTOOL_WARNING("DASTool", DASTOOL.out.das_bins_error.collect())
        ch_DBO_bin = DASTOOL.out.das_bins
        
        ch_checkm2_in = DASTOOL.out.das_bins
        
        multi_binner_log = BINNER_WARNING.out.log.mix(DASTOOL_WARNING.out.log)

        // new multi binner refinement tool
        MAGSCOT(ch_magscot_in,ch_magscot_folder)

    //*****************
    //If only one binning tool is selected, do not perform optimization.
    //*****************
    }else{
        //The input for checkm2 is the output folder of the selected binning tool.
        ch_checkm2_in = single_binner_result
    }
    


    //During quality control, if dastool encounters an exception, a log will be output, and no checkm2 task will be generated.
    CHECKM2 (ch_checkm2_in, ch_checkm2_db)
    CHECKM2_WARNING("CheckM2", CHECKM2.out.checkm2_error.collect())
    checkm2_log = CHECKM2_WARNING.out.log
    ch_DBO_qs = CHECKM2.out.quality_report

    QUAST_SAMPLE_BINS(ch_checkm2_in)

    //pick best bins result
    //bins floders channel [id, [DBO bins floder, PBO bins floders]]
    ch_best_bin = ch_DBO_bin.mix(ch_PBO_bin).groupTuple()

    //bins qs channel [id, [DBO bins qs list, PBO bins qs list]]
    ch_best_bin_qs = ch_DBO_qs.mix(ch_PBO_qs).groupTuple()


    ch_pcik_in = ch_best_bin.join(ch_best_bin_qs)

    PICKPBODBO(optimizeBins_method, ch_pcik_in) 
    bestBin = PICKPBODBO.out.bestBin

    //val(id),path(bins),path(quality_report)
    bestBin_report = PICKPBODBO.out.best

    binner_log = multi_binner_log.mix(checkm2_log)

    emit:
    contig_bowtie2
    bestBin
    bestBin_report
    binner_log

}

workflow BINNERCLEANUP {
    take:
    sample_number
    bestBin        //PICKPBODBO.out.bestBin
    bestBin_report  //PICKPBODBO.out.best

    main:

    // Mapping relationship between sampleID and bin filenames, using the channel merged from the final DBO / PBO.
    GETSAMPLEBINMAP(bestBin)
    sample_bin_map = GETSAMPLEBINMAP.out.mapping.collectFile(name: 'all_sample_bin_mapping.txt')

    //filter
    //filter_bins_input = ch_checkm2_in.join(CHECKM2.out.quality_report,by:0)
    FILTERBINS (bestBin_report)
    BIN_FILTER_WARNING("BinFilter", FILTERBINS.out.filtered_log.collect())
    filter_log = BIN_FILTER_WARNING.out.log

    //result 
    POSTBINNING (FILTERBINS.out.filter_bins.collect(), FILTERBINS.out.QS_quality_report.collect())
    binQS_report = POSTBINNING.out.rawbinQS
    all_filtered_bins = POSTBINNING.out.final_bins
    qs_quality_report = POSTBINNING.out.qs_quality_report

    // If ch_filtered_pass = 0, it means that all samples have no bins; output a log reminder.
    ch_filtered_pass = FILTERBINS.out.filter_bins.collect().count()
    NOBINSWARNING(ch_filtered_pass)
    nobins_log = NOBINSWARNING.out.log
    
    //Process input based on the number of samples.
    def ch_rename_bin = channel.empty()
    //single sample
    if (sample_number == 1){

        println "Input samplesheet has one sample."

        ch_rename_bin = all_filtered_bins

    }else{
        //drep
        DREP(all_filtered_bins)
        ch_rename_bin = DREP.out.genomes

        bestBin_report.collect().view()

        GALAH(all_filtered_bins,bestBin_report.map{ it[1] }.collect())
    }
    
    // Rename bins and split every 500 bins into a folder (to facilitate executing gtdb tasks).
    RENAMEBIN(ch_rename_bin)
    binInfo_report = RENAMEBIN.out.binInfo
    bins_folder = RENAMEBIN.out.folder
    final_genomes = RENAMEBIN.out.genomes
    bins_rename_map = RENAMEBIN.out.name_map
    bins_list = RENAMEBIN.out.list
    bins_count = RENAMEBIN.out.count
    bins_info = RENAMEBIN.out.bins_info

    binner_report = binQS_report.mix(binInfo_report, filter_log, nobins_log)

    emit:
    
    qs_quality_report
    bins_folder
    final_genomes
    bins_rename_map 
    bins_list
    bins_count
    bins_info
    binner_report
    sample_bin_map

}