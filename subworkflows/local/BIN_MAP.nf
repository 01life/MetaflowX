include { GETBINSAMPLE } from '../../modules/local/binassembly/bra_get_reassembly_bin_sample'

workflow BINMAP {
    take:
    clean_reads //ch_reads
    bin_QS_taxonomy //MERGEGTDB.out.bin_QS_taxonomy
    bin_count_abundance //06.BinsetProfile/061.BinAbundance/MetaFlowX_CoverM_bins_count_rename.xls
    bin_mean_abundance //6.BinsetProfile/061.BinAbundance/MetaFlowX_CoverM_bins_trimmed_mean_rename.xls
    bins_final_genomes //HQUniqueBins/folder*/*.fa //BINNER.out.final_genomes
    
    main:
    //=======================================================================================================//
    // ######## Step 1: Obtain refined & reassembled bin ID and sample ID, then merge reads channel. ########//
    //======================================================================================================//

        //select target bin and smaple to refine or reassembly
        GETBINSAMPLE(bin_QS_taxonomy, bin_count_abundance, bin_mean_abundance)
        binsample = GETBINSAMPLE.out.binsample

        //Construct ch_bin_fa, channel data structure: [bin, binfa].
        ch_bin_fa = bins_final_genomes.flatten().map { file ->
            def fileName = file.name.replaceFirst(/\.fa$/, '')
                [fileName, file]
            }

        //Construct ch_bin_sample, channel data structure: [bin, sample].
        ch_bin_sample = GETBINSAMPLE.out.binsample.splitText()
                .flatMap { line ->
                    def line_list = line.split("\t")
                    def bin = line_list[0].trim()
                    def contents = line_list[1].trim().split(',').collect { sample -> [ bin, sample] }
                }
        //Construct ch_bin_sample, channel data structure: [bin, sample, ref, nativesample].
        ch_bin_sample_ref = GETBINSAMPLE.out.binsample.splitText()
                .flatMap { line ->
                def line_list = line.split("\t")
                def bin = line_list[0].trim()
                def sample = line_list[1].trim()
                def ref = line_list[2].trim()
                def nativeSample = line_list[3].trim()
                return[[bin,sample,ref,nativeSample]]
                }
        // Reconstruct clean reads channel, channel data structure: [sampleID, sampleID, reads1, reads2] // To combine with ch_bin_sample channel by: 1.
        ch_clean_reads_formatted = clean_reads.map{ it -> [it[0],it[0],it[1][0],it[1][1] ]   }

        // Add reads information: combine ch_bin_sample with ch_clean_reads_formatted based on the second column. Original channel: [sampleID, binID, sampleID, reads1, reads2] reconstructed to: [binID, sampleID, reads1, reads2].
        ch_bin_sample_reads = ch_bin_sample.combine(ch_clean_reads_formatted, by:1).map { it -> [ it[1], it[0], it[3], it[4]]} 
        
        // Group by reads based on binID.
        ch_groupby_bin_reads = ch_bin_sample_reads.groupTuple().map{it -> [it[0],it[2],it[3]]}
        
        // Merge ch_bin_sample_ref & ch_groupby_bin_reads, channel data structure: [binID, sample_txt(,), ref, nativesample, [reads1, 2 list]].
        ch_bin_sample_ref_reads = ch_bin_sample_ref.join(ch_groupby_bin_reads)

        // cch_bin_sample_ref_reads & ch_bin_fa channel data structure: [binID, sample_txt(,), ref, nativesample, [reads1, 2 list], bin_fa].
        ch_bin_sample_ref_reads_binfa = ch_bin_sample_ref_reads.join(ch_bin_fa)


    emit:
    ch_bin_sample_ref_reads_binfa
    binsample


}