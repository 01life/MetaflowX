
//
// bin analysis
//

include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'

include { MERGEBINABUNTAXON } from '../../modules/local/binning/merge_bin_abun_taxon'

include { GENESET } from './GENESET'

include { BINNER } from './BINNER'
include { BINNERCLEANUP } from './BINNER'
include { BINTAXONOMY } from './BIN_TAXONOMY'
include { BINFUNCTION } from './BIN_FUNCTION'
include { BINABUNDANCE } from './BIN_ABUNDANCE'
include { BINMAP } from './BIN_MAP'
include { BINREFINE } from './BIN_REFINE'
include { BINREASSEMBLY } from './BIN_REASSEMBLY'


workflow BINNING {

    take:
    sample_number
    contigs           // channel: [ val(id), path(contigs) ]
    clean_reads       // channel: [ val(id), [ reads1, reads2 ] ]
    samplesheet
    prodigal_faa
    gene_info
    cdhit_clstr
    annotation
    vfdb_anno
    card_anno

    main:

    binning_report = Channel.empty()
    binning_bowtie2 = Channel.empty()

    bestBin =Channel.empty()
    bestBin_report = Channel.empty()
    bins_final_genomes = Channel.empty()
    
    binner_log = Channel.empty()
    // Input high quality bins (BINNER and BINNERCLEANUP module completed )
    if(params.HQ_unique_bins){
        bins_final_genomes = Channel.fromPath(params.HQ_unique_bins+"/*.fa").collect()
    // Input raw bins (BINNER module completed)
    }else if(params.rawbin_info){
        frawbin = file(params.rawbin_info, checkIfExists: true)
        // parse csv 
        bestBin_report = Channel.from (frawbin)
            .splitCsv (header:true, sep:',')
            .map { row ->  
                // check the expected columns
                def expectedColumns = ['id', 'rawbin_folder', 'quality_report']   
                def missingColumns = expectedColumns.findAll { !row.containsKey(it) }
                if (missingColumns.size() > 0) {
                    exit 1, "Configuration error: Missing required columns in '--rawbin_info' parameter, verify: ${missingColumns.join(', ')}"
                }
                // check empty entries in expectedColumns
                for (column in expectedColumns) {  
                    if (row[column] == null || row[column].trim() == '') {  
                        exit(1, "Empty value found in column: ${column}")  
                    }  
                }  
                // get information
                return [row.id, row.rawbin_folder, row.quality_report]  
            }  
        bestBin = bestBin_report.map{ it -> [it[0], it[1]]}
    }else{
        //Bining
        BINNER(contigs, clean_reads, prodigal_faa)
        binning_bowtie2 = binning_bowtie2.mix(BINNER.out.contig_bowtie2)
        binner_log = BINNER.out.binner_log
        bestBin = BINNER.out.bestBin
        bestBin_report = BINNER.out.bestBin_report
    }

    binner_report = Channel.empty()
    bins_list = Channel.empty()
    bins_folder = Channel.empty()
    bins_count = Channel.empty()
    qs_quality_report = Channel.empty()
    bins_info = Channel.empty() 
    bins_rename_map = Channel.empty()
    sample_bin_map = Channel.empty()

    def require_hq_bins = params.bin_taxonomy || params.bin_function || params.bin_abundance_calculation
    if(require_hq_bins && !params.HQ_unique_bins){
        // Get HQ unique bins
        BINNERCLEANUP(sample_number, bestBin, bestBin_report)
        binner_report = BINNERCLEANUP.out.binner_report
        bins_list = BINNERCLEANUP.out.bins_list
        bins_folder = BINNERCLEANUP.out.bins_folder
        bins_count = BINNERCLEANUP.out.bins_count
        qs_quality_report = BINNERCLEANUP.out.qs_quality_report
        bins_info = BINNERCLEANUP.out.bins_info 
        bins_rename_map = BINNERCLEANUP.out.bins_rename_map
        bins_final_genomes = BINNERCLEANUP.out.final_genomes
        sample_bin_map = BINNERCLEANUP.out.sample_bin_map
    }

    //bin taxonomy
    bin_QS_taxonomy = Channel.empty()
    gtdb_result = Channel.empty()
    taxonomy_report= Channel.empty()
    if(params.bin_taxonomy){
        BINTAXONOMY(bins_folder, bins_rename_map, qs_quality_report, sample_bin_map)
        gtdb_summary = BINTAXONOMY.out.gtdb_summary.map{
            def summary = file(it)
            return "${summary.parent}/${summary.name}";
        }
        bin_QS_taxonomy = BINTAXONOMY.out.bin_QS_taxonomy
        gtdb_result = BINTAXONOMY.out.gtdb_summary
        taxonomy_report = BINTAXONOMY.out.gtdb_report
    }

    //bin gene function
    function_report = Channel.empty()
    if(params.bin_function){
        BINFUNCTION(gene_info, cdhit_clstr, annotation, bins_final_genomes, vfdb_anno, card_anno)
        function_report = BINFUNCTION.out.binFunction_geneID_report
    }

    //Bin Abundance Calculation.
    depth_list = Channel.empty()
    bins_rel_abun = Channel.empty()
    bins_mean_abun = Channel.empty()
    abundance_report = Channel.empty()
    if(params.bin_abundance_calculation){
        ch_method4coverm = params.method4coverm
        if(params.bin_reassembly){
            if(!params.method4coverm.contains("relative_abundance")){
                ch_method4coverm = ch_method4coverm.concat(",relative_abundance")
            }
            if(!params.method4coverm.contains("trimmed_mean")){
                ch_method4coverm = ch_method4coverm.concat(",trimmed_mean")
            }
            if(!params.method4coverm.contains("count")){
                ch_method4coverm = ch_method4coverm.concat(",count")
            }
        }
        BINABUNDANCE(sample_number, clean_reads, bins_final_genomes, ch_method4coverm, samplesheet)
        abundance_report = BINABUNDANCE.out.binAundance_report
        binning_bowtie2 = binning_bowtie2.mix(BINABUNDANCE.out.bin_bowtie2)
        depth_list = BINABUNDANCE.out.depth_list
        bins_rel_abun = BINABUNDANCE.out.totalRelativeAbun
        bins_count_abun = BINABUNDANCE.out.totalCountAbun
        bins_mean_abun = BINABUNDANCE.out.totalMeanAbun
    }
    
    // Bin Taxonomy Classification.
    ch_bins_taxon = Channel.empty()
    MERGEBINABUNTAXON(bin_QS_taxonomy, bins_rel_abun)
    ch_bins_taxon = MERGEBINABUNTAXON.out.taxon
    
    //Bin Refinement and Reassembly.
    ch_filter_bin  = Channel.empty()
    ch_filter_bin_info = Channel.empty()
    ch_filter_bin_mash_fq = Channel.empty()

    //Only supports PE data.
    rebin_report = Channel.empty()
    if( !params.single_end && (params.bin_refine || params.bin_reassembly) ){
        
        //bin map
        BINMAP(clean_reads, bin_QS_taxonomy, bins_count_abun, bins_mean_abun, bins_final_genomes)
        ch_pick_info = BINMAP.out.ch_bin_sample_ref_reads_binfa
        ch_binsample = BINMAP.out.binsample

        // bin refine
        if(params.bin_refine){
            BINREFINE(clean_reads, contig, ch_filter_bin, ch_filter_bin_info, ch_filter_bin_mash_fq)
            ch_filter_bin_mash_fq = BINREFINE.out.after_refine_bin_mash_fq
        }

        //bin reassembly
        if(params.bin_reassembly){
            BINREASSEMBLY(ch_pick_info, ch_binsample, clean_reads, bin_QS_taxonomy)
            rebin_report = BINREASSEMBLY.out.rebin_report
        }
    }
    
    binning_report = binning_report.mix(binner_log, binner_report, taxonomy_report, function_report, abundance_report, rebin_report).collect()

    emit:
    bins_list
    bins_folder
    bins_count
    qs_quality_report
    bins_info
    bins_rename_map
    gtdb_result
    depth_list
    bins_rel_abun
    binning_report
    binning_bowtie2
   
}