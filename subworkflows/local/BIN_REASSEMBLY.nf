
include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'

include { EXTRACTREADSMASHSAMPLE } from '../../modules/local/binassembly/bra_extract_bin_mapReads_get_bestMash_sample'
include { REBINASSEMBLY } from '../../modules/local/binassembly/bra_reBinAssembly'
include { MEGAHIT as BRAMEGAHIT } from '../../modules/local/binassembly/bra_megahit'
include { BRACONTIGFILTER } from '../../modules/local/binassembly/bra_contig_filter'
include { REBINSTAT } from '../../modules/local/binassembly/bra_rebin_stat'
include { CHECKM2 as ASS_CHECKM2 } from '../../modules/local/binning/checkm2'
include { FILTERBINS as ASS_FILTERBINS} from '../../modules/local/binning/filter_bins'
include { PIPELINEWARNING as ASS_CHECKM2_WARNING } from '../../modules/local/common/pipeline_warning'
include { OUTHQBIN } from '../../modules/local/binassembly/bra_pick_reassembly_hq_bin'
include { BRADEPTH } from '../../modules/local/binassembly/bra_depth'
include { MULTICOV } from '../../modules/local/binassembly/bra_merge_multi_depth'
include { BRAMETABAT2 } from '../../modules/local/binassembly/bra_metabat2'
include { BRASEMIBIN2 } from '../../modules/local/binassembly/bra_semibin2'
include { BRADASTOOL } from '../../modules/local/binassembly/bra_das_tool'
include { PIPELINEWARNING as BRADASTOOL_WARNING } from '../../modules/local/common/pipeline_warning'
include { CHECKM2 as R2_CHECKM2 } from '../../modules/local/binning/checkm2'
include { FILTERBINS as R2_FILTERBINS} from '../../modules/local/binning/filter_bins'
include { PIPELINEWARNING as R2_CHECKM2_WARNING } from '../../modules/local/common/pipeline_warning'
include { PICKREBIN } from '../../modules/local/binassembly/bra_pick_rebin_HQ_bin'
include { MERGEREBINCHECKM2 } from '../../modules/local/binassembly/bra_merge_rebin_checkm2'
include { MERGEREBINUNIMPROVEBIN } from '../../modules/local/binassembly/bra_merge_rebin_unimprove'
include { DEEPURIFYCLEAN } from '../../modules/local/binassembly/bra_deepurify_clean'
include { DEEPURIFYCLEANRENAME } from '../../modules/local/binassembly/bra_deepurify_clean'
include { PICKREREFINE } from '../../modules/local/binassembly/bra_pick_re-refine_HQ_bin'
include { REBINVERIFY } from '../../modules/local/binassembly/bra_rebin_verify_noRefine'


workflow BINREASSEMBLY {
    take:
    ch_bin_sample_ref_reads //[binis,sample_txt(,),ref,nativesample,[reads1,2 list],bin_fa]
    binsample
    clean_reads
    ch_bin_QS_taxonomy


    main:

    /*
    * Verify the essential parameters for running this module
    */
    binner_essential_db = [params.checkm2_db]

    if(!checkEssentialParams(binner_essential_db)) { exit 1, "The required parameter to execute the Bin Reassembly module is: checkm2_db" }

    ch_checkm2_db = params2Channel(params.checkm2_db)


    // Extract map reads.
    EXTRACTREADSMASHSAMPLE(ch_bin_sample_ref_reads)
    // Flatten the sample IDs that match the target bin, channel data structure: [binid, sampleID, binid-sampleID] // One bin ID can correspond to multiple samples.
    ch_bin_sample = EXTRACTREADSMASHSAMPLE.out.binsample.collectFile(name: 'target_bin_best_mash_sample.txt').splitText()
            .map { line ->
                def line_list = line.split("\t")
                def bin = line_list[0].trim()
                def sample = line_list[1].trim()
                return [ bin,sample,bin + "-" + sample ]
                 }
            

    // re-assembly
    ch_reassembly_input = ch_bin_sample_ref_reads.map{ it -> [it[0], it[6]] }.join(EXTRACTREADSMASHSAMPLE.out.map)
    REBINASSEMBLY(ch_reassembly_input)

    // Switch to MEGAHIT if sapdes assembly fails.
    spades_contigs = REBINASSEMBLY.out.rebin
    ch_spades_failed = ch_reassembly_input.join(spades_contigs, remainder: true).filter{ it[2] == null }.map{ it -> [it[0], it[1]]}
    
    BRAMEGAHIT(ch_spades_failed)
    megahit_contigs = BRAMEGAHIT.out.contigs

    contigs = spades_contigs.concat(megahit_contigs)

    // Filter length & rename.
    BRACONTIGFILTER(contigs)

    // Contig statistics.
    REBINSTAT("reassembly", BRACONTIGFILTER.out.filtercontigs.collect())

    // First check quality control (checkm2) after assembly.
    // ASSCHECKM2(outdir, REBINSTAT.out.reassembly_bins)

    ch_checkm2_in = REBINSTAT.out.reassembly_bins
    ASS_CHECKM2(ch_checkm2_in, ch_checkm2_db)
    ASS_CHECKM2_WARNING("BRA_ASS_CheckM2", ASS_CHECKM2.out.checkm2_error.collect())

    //filter
    filter_bins_input = ch_checkm2_in.join(ASS_CHECKM2.out.quality_report,by:0)
    ASS_FILTERBINS (filter_bins_input)


    // First round of high-quality bin release.
    OUTHQBIN('R1_After_Reassembly_HQ_Bins',BRACONTIGFILTER.out.filtercontigs.collect(),ASS_FILTERBINS.out.QS_quality_report)

    // // Final publication of bin information for the first round of high-quality bin release.
    // // OUTHQBIN.out.improved_info

    // Generate first round of low-quality bin fa channel, channel data structure: [binid, binfa] // Split binid by "_" and take the first element.
    ch_LQ_bin = OUTHQBIN.out.lqbin.collect()
                        .flatten()
                        .map{ binfa ->
                            [ binfa.baseName.split('_reassembly_contigs_')[0] , binfa ]
                        }
    
    // Reconstruct clean reads channel, channel data structure: [sampleID, sampleID, reads1, reads2] // To combine with ch_bin_sample channel by: 1.
    ch_clean_reads_formatted = clean_reads.map{ it -> [it[0],it[0],it[1][0],it[1][1]]   }

    // Add reads information: combine ch_bin_sample with ch_clean_reads_formatted based on the second column. Original channel: [sampleID, binID, binID-sampleID, sampleID, reads1, reads2] reconstructed to: [binID, sampleID, binID-sampleID, reads1, reads2].
    ch_bin_sample_reads = ch_bin_sample.combine(ch_clean_reads_formatted, by:1).map { it -> [ it[1], it[3], it[2], it[4], it[5] ]} 

    // Add binfa information: combine ch_bin_sample_reads with ch_LQ_bin based on the first column, channel data structure: [binID, sampleID, binID-sampleID, reads1, reads2, binfa].
    ch_bin_sample_reads_binfa = ch_bin_sample_reads.combine(ch_LQ_bin,by: 0)

    // Calculate depth for each binID-sampleID, output data structure: out = [binID, sampleID, binID-sampleID, bam/depth].
    BRADEPTH(ch_bin_sample_reads_binfa)

    // Group depth files by binID.
    BRADEPTH.out.bin_depth.groupTuple().map{it -> [it[0], it[3] ] }

    // Merge multiple sample depth files for bins into one depth matrix.
    MULTICOV(BRADEPTH.out.bin_depth.groupTuple().map{it -> [it[0], it[3] ]})


    // LQ bin fa join LQ bin  multi depth
    ch_LQ_bin_multiCov = ch_LQ_bin.join(MULTICOV.out.mulitCov)

    // Multi-cover with matebat2 for binning, output data structure: out = [binID, metabat.contigs2bin.tsv] // One bin produces one matebat2 result.
    BRAMETABAT2(ch_LQ_bin_multiCov)

    // Combine depth bam files for ch_bin_sample_reads_binfa based on the third column binID-sampleID, original channel： [0:binID-sampleID, 1:binID, 2:sampleID, 3:reads1, 4:reads2, 5:binfa, 6:binID, 7:sampleID, 8:bam ] econstructed to ：[0:binID, 1:sampleID, 2:binID-sampleID, 3:reads1, 4:reads2, 5:binsfa, 6:bam]
    ch_bin_sample_reads_binfa_bam = ch_bin_sample_reads_binfa.combine(BRADEPTH.out.sorted_bam, by:2).map{ it -> [ it[1], it[2], it[0], it[3], it[4] , it[5], it[8] ]}

    // For each binID-sampleID, re-bin with semibin2.
    BRASEMIBIN2(ch_bin_sample_reads_binfa_bam.map{ it -> [it[0], it[1], it[5], it[6] ]})


    //Group results of semibin2 for each binID-sampleID, channel data structure: [binid, [semebin2.tsv1, ...semebin2.tsvN]].
    ch_semibin_csv = BRASEMIBIN2.out.tsv.groupTuple().map{it -> [it[0], it[2] ]}

    //Merge results of metabat2 & semibin2, channel data structure: [binid, [metabat2.tsv, semebin2.tsv1, ...semebin2.tsvN], binfa].
    multi_binner_result = BRAMETABAT2.out.tsv.combine(ch_semibin_csv, by: 0)
                 .map { item ->
                     def (binID, metabat2, semibinlist) = item
                     return [binID, [metabat2] + semibinlist]
                 }.combine(ch_LQ_bin, by:0)
    // Select optimal results from metabat2 & semibin2.
    BRADASTOOL(multi_binner_result)


    // After re-binning, perform quality control (checkm2) on the selected results.
    ch_checkm2_in_r2  = BRADASTOOL.out.das_bins
    R2_CHECKM2(ch_checkm2_in_r2, ch_checkm2_db)
    R2_CHECKM2_WARNING("BRA_R2_CheckM2", R2_CHECKM2.out.checkm2_error.collect())

    //filter
    // filter_bins_input = ch_checkm2_in_r2.join(R2_CHECKM2.out.quality_report,by:0)
    // R2_FILTERBINS (filter_bins_input)


    BRADASTOOL.out.das_bins
        .filter { it -> it != null && !it.isEmpty() }
        .set { valid_das_bins }

    // valid_das_bins_error = BRADASTOOL.out.das_bins_error.collect()
    //     .filter { it -> it != null && !it.isEmpty() }
    //     .ifEmpty { Channel.empty() }
    
    // valid_das_bins_error.collect().view()
    // BRADASTOOL_WARNING("BRA_DASTool", valid_das_bins_error.collect())
    BRADASTOOL.out.das_bins_error.collect().view()
    BRADASTOOL_WARNING("BRA_DASTool", BRADASTOOL.out.das_bins_error.collect())

    
    //Select the highest quality bin and release it.
    ch_pcik_rebin = R2_CHECKM2.out.quality_report.join(valid_das_bins).combine(ch_bin_QS_taxonomy)
    // ch_pcik_rebin.view()

    // Second round of high-quality bin release.
    PICKREBIN(ch_pcik_rebin)
    // ReAss_ReBin_HQ  ReAss_ReBin_Max_QS  ReAss_ReBin_Unimprove

    rebinfa_branched = PICKREBIN.out.rebinfa
        .branch {
            Max_QS: it[1].name.contains('Max_QS')
                return it
            Unimprove: it[1].name.contains('Unimprove')
                return it
                
        }

    rebinfa_branched_unimprove = rebinfa_branched.Unimprove.map{ it -> it[1] }.collect()

    MERGEREBINCHECKM2(PICKREBIN.out.rebin_org_Checkm2.collect(),PICKREBIN.out.improved_info.collect())

    // Final publication of bin information for the second round of high-quality bin release.
    // MERGEREBINCHECKM2.out.improve_info

    MERGEREBINUNIMPROVEBIN(rebinfa_branched_unimprove)

    // Re-bin unselected bins for deep purification.

    DEEPURIFYCLEAN(rebinfa_branched.Max_QS, ch_checkm2_db)

    DEEPURIFYCLEANRENAME(DEEPURIFYCLEAN.out.deepurifyDir)
    

    // // // Third round of high-quality bin release.
    ch_merge_rebin_checkm2 = PICKREBIN.out.rebinCheckm2
        .collectFile(name: 'merged_ReAss_ReBin_checkm2.txt') // newLine: true


    PICKREREFINE(DEEPURIFYCLEANRENAME.out.renameDeepurify.map{it -> it[1]}.collect(), ch_merge_rebin_checkm2)

    // Final publication of bin information for the third round of high-quality bin release.
    // PICKREREFINE.out.improved_info

    // OUTHQBIN.out.improved_info.view()
    // ch_merge_improve = OUTHQBIN.out.improved_info.mix(MERGEREBINCHECKM2.out.improved_info, PICKREREFINE.out.improved_info).collectFile(name: 'merged_HQ_improve.txt') // newLine: true

    // REBINVERIFY(binsample,ch_bin_QS_taxonomy,OUTHQBIN.out.reass_checkm2,MERGEREBINCHECKM2.out.rebining_checkm2,PICKREREFINE.out.refine_checkm2, ch_merge_improve)

    REBINVERIFY(binsample,ch_bin_QS_taxonomy,OUTHQBIN.out.reass_checkm2,MERGEREBINCHECKM2.out.rebining_checkm2,PICKREREFINE.out.refine_checkm2, OUTHQBIN.out.improved_info, MERGEREBINCHECKM2.out.improved_info, PICKREREFINE.out.improved_info)
    rebin_report = REBINVERIFY.out.bin_reassembly_stat_info

    emit:
    rebin_report

}