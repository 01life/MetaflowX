/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetassembly.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist

def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters

if (!(params.mode in [0,1,2,3,4,5])) { exit 1, "The parameter mode is invalid, supported values are 0, 1, 2, 3, 4, 5 !" }

if (params.outdir) { ch_output = new File(params.outdir).getAbsolutePath()  } else { exit 1, "Output directory not specified !" }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { params2Channel } from '../modules/local/common/utils'
include { checkEssentialParams } from '../modules/local/common/utils'
include { PRODIGAL } from '../modules/local/geneset/prodigal'
include { CONTIGFILTER } from '../modules/local/common/contig_filter'
include { CONTIGSTAT } from '../modules/local/assembly/contig_stat'
include { MERGEBINABUNTAXON } from '../modules/local/binning/merge_bin_abun_taxon'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/INPUT_CHECK'
include { QC } from '../subworkflows/local/QC'
include { ASSEMBLY } from '../subworkflows/local/ASSEMBLY'
include { RAPID_TAXONOMIC_PROFILING } from '../subworkflows/local/RAPID_TAXONOMIC_PROFILING'
include { GENESET } from '../subworkflows/local/GENESET'
include { BINNER } from '../subworkflows/local/BINNER'
include { BINTAXONOMY } from '../subworkflows/local/BIN_TAXONOMY'
include { BINFUNCTION } from '../subworkflows/local/BIN_FUNCTION'
include { BINABUNDANCE } from '../subworkflows/local/BIN_ABUNDANCE'
include { BINMAP } from '../subworkflows/local/BIN_MAP'
include { BINREFINE } from '../subworkflows/local/BIN_REFINE'
include { BINREASSEMBLY } from '../subworkflows/local/BIN_REASSEMBLY'
include { POLISH } from '../subworkflows/local/POLISH'


////////////////////////////////////////////////////
/* -----   Define channel for SUBWORKFLOW -----   */
////////////////////////////////////////////////////      

ch_raw_reads    = Channel.empty()
ch_clean_reads  = Channel.empty()

ch_prodigal_faa = Channel.empty()
ch_cdhit_clstr  = Channel.empty()
ch_gene_info    = Channel.empty()
ch_annotation   = Channel.empty()
ch_vfdb_anno    = Channel.empty()
ch_card_anno    = Channel.empty()

//html report input
ch_report_input = Channel.empty()
ch_fastp_json = Channel.empty()
ch_bowtie2_log = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary


workflow METASSEMBLY {
    
    def ch_input = null 

    if (params.input) { ch_input = file(params.input) } else { exit 1, "Input samplesheet not specified !" }
    
    // total sample number
    def sample_number = ch_input.splitCsv().size() - 1

    //Input Check and Data Parsing
    
    //Determine if parsing is for raw_reads/clean_reads based on column names.
    def columns = ch_input.text.tokenize('\n').first()
    def raw_reads_flag = false

    //Perform QC and input raw_reads when executing the full process.
    if( params.mode==1 || (params.mode==0 && !params.skip_qc) ){
        // println("raw_reads")
        if(columns.contains("raw")) {
            INPUT_CHECK("raw", ch_input)
        }else{ exit 1, "Invalid input samplesheet: expects column raw_reads1 and raw_reads2 or raw_se !"}
    //Modes 2 and 3 can accept raw reads.
    }else if( params.mode in [2, 3] ){
        if(columns.contains("clean")) {
            INPUT_CHECK("clean", ch_input)
        }else if(columns.contains("raw")){
            INPUT_CHECK("raw", ch_input)
            raw_reads_flag = true
        }else{ exit 1, "Invalid input samplesheet: expects column raw_reads1,raw_reads2,raw_se or clean_reads1,clean_read2,clean_se !"}
    //For inputting clean reads.
    }else{
        // println("clean_reads")
        if(columns.contains("clean")) { 
            INPUT_CHECK("clean", ch_input)
        }else{ exit 1, "Invalid input samplesheet: expects column clean_reads !"}
    }
    
    ch_raw_reads = INPUT_CHECK.out.raw_reads

    ch_clean_reads = INPUT_CHECK.out.clean_reads        
    
    CONTIGFILTER(INPUT_CHECK.out.contig)
    ch_contig = CONTIGFILTER.out.contigs

    //Submodule Execution.
    switch(params.mode){
        //QC
        case 1 : {
            QC(sample_number, ch_raw_reads)
            ch_report_input = ch_report_input.mix(QC.out.qc_report)
            ch_clean_reads = QC.out.clean_reads
            ch_fastp_json = QC.out.fastp_json
            break;
        }
        //Assembly.
        case 2 : {
            
            //Input raw reads, perform QC first.
            if(raw_reads_flag){
                QC(sample_number, ch_raw_reads)
                ch_report_input = ch_report_input.mix(QC.out.qc_report)
                ch_clean_reads = QC.out.clean_reads
                ch_fastp_json = QC.out.fastp_json
            }

            ASSEMBLY(sample_number, ch_clean_reads)
            ch_report_input = ch_report_input.mix(ASSEMBLY.out.assembly_report)
            break;
        }
        //reads marker
        case 3 : {
            //Input raw reads, perform QC first.
            if(raw_reads_flag){
                QC(sample_number, ch_raw_reads)
                ch_report_input = ch_report_input.mix(QC.out.qc_report)
                ch_clean_reads = QC.out.clean_reads
                ch_fastp_json = QC.out.fastp_json
            }

            RAPID_TAXONOMIC_PROFILING(sample_number, ch_clean_reads, ch_input)
            ch_report_input = ch_report_input.mix(RAPID_TAXONOMIC_PROFILING.out.marker_report)
            ch_marker_profile = RAPID_TAXONOMIC_PROFILING.out.ch_profile

            break;
        }
        //Gene Set Prediction and Annotation.
        case 4 : {
            GENESET(sample_number, ch_contig, ch_clean_reads, ch_input)
            ch_report_input = ch_report_input.mix(GENESET.out.geneset_report)
            ch_bowtie2_log = ch_bowtie2_log.mix(GENESET.out.geneset_bowtie2)
            ch_geneset_profile = GENESET.out.ch_profile
            break;
        }
        //Binning and Bins Annotation.
        case 5 : {

            CONTIGSTAT(sample_number, sample_number, ch_contig.map{ _,file -> file}.collect())
            ch_contig_info = CONTIGSTAT.out.contig_info

            //Bining
            PRODIGAL (ch_contig)
            BINNER(sample_number, ch_contig, ch_clean_reads, PRODIGAL.out.faa)
            ch_report_input = ch_report_input.mix(BINNER.out.binner_report)
            ch_bowtie2_log = ch_bowtie2_log.mix(BINNER.out.contig_bowtie2)
            ch_bins_list = BINNER.out.bins_list
            ch_bins_folder = BINNER.out.bins_folder
            ch_bins_count = BINNER.out.bins_count
            ch_qs_quality_report = BINNER.out.qs_quality_report
            ch_bins_info = BINNER.out.bins_info 
            ch_bins_rename_map = BINNER.out.bins_rename_map
            ch_bins_final_genomes = BINNER.out.final_genomes
            ch_sample_bin_map = BINNER.out.sample_bin_map

            //bin taxonomy
            ch_gtdb_summary = "N"
            ch_bin_QS_taxonomy = Channel.empty()
            if(params.bin_taxonomy){
                BINTAXONOMY(ch_bins_folder, ch_bins_rename_map, ch_qs_quality_report, ch_sample_bin_map)
                ch_gtdb_summary = BINTAXONOMY.out.gtdb_summary.map{
                    def summary = file(it)
                    return "${summary.parent}/${summary.name}";
                }
                ch_report_input = ch_report_input.mix(BINTAXONOMY.out.gtdb_report) 
                ch_gtdb_result = BINTAXONOMY.out.gtdb_summary
                ch_bin_QS_taxonomy = BINTAXONOMY.out.bin_QS_taxonomy
            }

            //bin gene function
            if(params.bin_function){
                
                //Check Input Parameters.
                extraParamsList = [ params.gene_info, params.cdhit_clstr, params.emapper_annotation ]
                extra_flag = checkEssentialParams(extraParamsList)
                
                //When the required configuration files (three) exist, only execute prodigal; otherwise, execute the GENESET submodule.
                if(extra_flag){
                    ch_gene_info = params2Channel(params.gene_info)
                    ch_cdhit_clstr = params2Channel(params.cdhit_clstr)
                    ch_annotation = params2Channel(params.emapper_annotation)
                    ch_vfdb_anno = params2Channel(params.vfdb_annotation)
                    ch_card_anno = params2Channel(params.card_annotation)
                }else{
                    GENESET(sample_number, ch_contig, ch_clean_reads, ch_input)
                    ch_prodigal_faa = GENESET.out.prodigal_faa
                    ch_cdhit_clstr = GENESET.out.clstr
                    ch_gene_info = GENESET.out.gene_info
                    ch_annotation = GENESET.out.annotation
                    ch_vfdb_anno = GENESET.out.vfdb_anno
                    ch_card_anno = GENESET.out.card_anno
                    
                    ch_report_input = ch_report_input.mix(GENESET.out.geneset_report)
                    ch_bowtie2_log = ch_bowtie2_log.mix(GENESET.out.geneset_bowtie2)
                    ch_geneset_profile = GENESET.out.ch_profile
                
                }

                BINFUNCTION(ch_gene_info, ch_cdhit_clstr, ch_annotation, ch_bins_final_genomes,ch_vfdb_anno,ch_card_anno)
                ch_report_input = ch_report_input.mix(BINFUNCTION.out.binFunction_geneID_report)

            }

            //Bin Abundance Calculation.
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
                BINABUNDANCE(sample_number, ch_clean_reads, ch_bins_final_genomes, ch_bins_folder, ch_method4coverm, ch_input)
                ch_report_input = ch_report_input.mix(BINABUNDANCE.out.binAundance_report) 
                ch_bowtie2_log = ch_bowtie2_log.mix(BINABUNDANCE.out.bin_bowtie2)
                ch_depth_list = BINABUNDANCE.out.depth_list
                ch_bins_rel_abun = BINABUNDANCE.out.totalRelativeAbun
                ch_bins_count_abun = BINABUNDANCE.out.totalCountAbun
                ch_bins_mean_abun = BINABUNDANCE.out.totalMeanAbun
            }

            // Bin Taxonomy Classification.
            ch_bins_taxon = Channel.empty()
            MERGEBINABUNTAXON(ch_bin_QS_taxonomy, ch_bins_rel_abun)
            ch_bins_taxon = MERGEBINABUNTAXON.out.taxon


            //Bin Refinement and Reassembly.
            ch_filter_bin  = Channel.empty()
            ch_filter_bin_info = Channel.empty()
            ch_filter_bin_mash_fq = Channel.empty()

            //Only supports PE data.
            if( !params.single_end && (params.bin_refine || params.bin_reassembly) ){
                
                //bin map
                BINMAP(ch_clean_reads, ch_bin_QS_taxonomy, ch_bins_count_abun, ch_bins_mean_abun, ch_bins_final_genomes)
                ch_pick_info = BINMAP.out.ch_bin_sample_ref_reads_binfa
                ch_binsample = BINMAP.out.binsample

                // bin refine
                if(params.bin_refine){
                    BINREFINE(ch_clean_reads, ch_contig, ch_filter_bin, ch_filter_bin_info, ch_filter_bin_mash_fq)
                    ch_filter_bin_mash_fq = BINREFINE.out.after_refine_bin_mash_fq
                }

                //bin reassembly
                if(params.bin_reassembly){
                    BINREASSEMBLY(ch_pick_info, ch_binsample, ch_clean_reads, ch_bin_QS_taxonomy)
                    ch_report_input = ch_report_input.mix(BINREASSEMBLY.out.rebin_report)
                }

            }

            break;

        }
        //Default executes the entire process.
        default: {
            
            if(!params.skip_qc){
                QC(sample_number, ch_raw_reads)
                ch_clean_reads = QC.out.clean_reads
                ch_report_input = ch_report_input.mix(QC.out.qc_report)
                ch_fastp_json = QC.out.fastp_json
            }

            ASSEMBLY(sample_number, ch_clean_reads)
            ch_contig = ASSEMBLY.out.contigs
            ch_report_input = ch_report_input.mix(ASSEMBLY.out.assembly_report)
            ch_contig_info = ASSEMBLY.out.contig_info

            if(!params.skip_marker){
                RAPID_TAXONOMIC_PROFILING(sample_number, ch_clean_reads, ch_input)
                ch_report_input = ch_report_input.mix(RAPID_TAXONOMIC_PROFILING.out.marker_report)
                ch_marker_profile = RAPID_TAXONOMIC_PROFILING.out.ch_profile
            }
            
            GENESET(sample_number, ch_contig, ch_clean_reads, ch_input)
            ch_prodigal_faa = GENESET.out.prodigal_faa
            ch_cdhit_clstr = GENESET.out.clstr
            ch_gene_info = GENESET.out.gene_info
            ch_annotation = GENESET.out.annotation
            ch_report_input = ch_report_input.mix(GENESET.out.geneset_report)
            ch_bowtie2_log = ch_bowtie2_log.mix(GENESET.out.geneset_bowtie2)
            ch_geneset_profile = GENESET.out.ch_profile
            ch_vfdb_anno = GENESET.out.vfdb_anno
            ch_card_anno = GENESET.out.card_anno
            
            if(!params.skip_binning){
                //Bining
                BINNER(sample_number, ch_contig, ch_clean_reads, ch_prodigal_faa)
                ch_report_input = ch_report_input.mix(BINNER.out.binner_report)
                ch_bowtie2_log = ch_bowtie2_log.mix(BINNER.out.contig_bowtie2)
                ch_bins_list = BINNER.out.bins_list
                ch_bins_folder = BINNER.out.bins_folder
                ch_bins_count = BINNER.out.bins_count
                ch_qs_quality_report = BINNER.out.qs_quality_report
                ch_bins_info = BINNER.out.bins_info 
                ch_bins_rename_map = BINNER.out.bins_rename_map
                ch_bins_final_genomes = BINNER.out.final_genomes
                ch_sample_bin_map = BINNER.out.sample_bin_map

                //bin taxonomy
                ch_gtdb_summary = "N"
                ch_bin_QS_taxonomy = Channel.empty()
                if(params.bin_taxonomy){
                    BINTAXONOMY(ch_bins_folder, ch_bins_rename_map, ch_qs_quality_report, ch_sample_bin_map)
                    ch_gtdb_summary = BINTAXONOMY.out.gtdb_summary.map{
                        def summary = file(it)
                        return "${summary.parent}/${summary.name}";
                    }
                    ch_report_input = ch_report_input.mix(BINTAXONOMY.out.gtdb_report) 
                    ch_gtdb_result = BINTAXONOMY.out.gtdb_summary
                    ch_bin_QS_taxonomy = BINTAXONOMY.out.bin_QS_taxonomy
                }

                //bin gene function
                if(params.bin_function){
                    BINFUNCTION(ch_gene_info, ch_cdhit_clstr, ch_annotation, ch_bins_final_genomes, ch_vfdb_anno, ch_card_anno)
                    ch_report_input = ch_report_input.mix(BINFUNCTION.out.binFunction_geneID_report)
                }

                //Bin Abundance Calculation.
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
                    BINABUNDANCE(sample_number, ch_clean_reads, ch_bins_final_genomes, ch_bins_folder, ch_method4coverm, ch_input)
                    ch_report_input = ch_report_input.mix(BINABUNDANCE.out.binAundance_report) 
                    ch_bowtie2_log = ch_bowtie2_log.mix(BINABUNDANCE.out.bin_bowtie2)
                    ch_depth_list = BINABUNDANCE.out.depth_list
                    ch_bins_rel_abun = BINABUNDANCE.out.totalRelativeAbun
                    ch_bins_count_abun = BINABUNDANCE.out.totalCountAbun
                    ch_bins_mean_abun = BINABUNDANCE.out.totalMeanAbun
                }
                
                // Bin Taxonomy Classification.
                ch_bins_taxon = Channel.empty()
                MERGEBINABUNTAXON(ch_bin_QS_taxonomy, ch_bins_rel_abun)
                ch_bins_taxon = MERGEBINABUNTAXON.out.taxon
                
                //Bin Refinement and Reassembly.
                ch_filter_bin  = Channel.empty()
                ch_filter_bin_info = Channel.empty()
                ch_filter_bin_mash_fq = Channel.empty()

                //Only supports PE data.
                if( !params.single_end && (params.bin_refine || params.bin_reassembly) ){
                    
                    //bin map
                    BINMAP(ch_clean_reads, ch_bin_QS_taxonomy, ch_bins_count_abun, ch_bins_mean_abun, ch_bins_final_genomes)
                    ch_pick_info = BINMAP.out.ch_bin_sample_ref_reads_binfa
                    ch_binsample = BINMAP.out.binsample

                    // bin refine
                    if(params.bin_refine){
                        BINREFINE(ch_clean_reads, ch_contig, ch_filter_bin, ch_filter_bin_info, ch_filter_bin_mash_fq)
                        ch_filter_bin_mash_fq = BINREFINE.out.after_refine_bin_mash_fq
                    }

                    //bin reassembly
                    if(params.bin_reassembly){
                        BINREASSEMBLY(ch_pick_info, ch_binsample, ch_clean_reads, ch_bin_QS_taxonomy)
                        ch_report_input = ch_report_input.mix(BINREASSEMBLY.out.rebin_report)
                    }

                }

            }
            
            break;
        
        }

    }

    // Visualization Report.
    POLISH(ch_report_input.collect(), ch_fastp_json, ch_bowtie2_log.collect(), ch_clean_reads, raw_reads_flag)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
