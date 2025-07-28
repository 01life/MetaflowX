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
include { CONTIGFILTER } from '../modules/local/common/contig_filter'
include { CONTIGSTAT } from '../modules/local/assembly/contig_stat'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/INPUT_CHECK'
include { QC } from '../subworkflows/local/QC'
include { ASSEMBLY } from '../subworkflows/local/ASSEMBLY'
include { RAPID_TAXONOMIC_PROFILING } from '../subworkflows/local/RAPID_TAXONOMIC_PROFILING'
include { GENESET } from '../subworkflows/local/GENESET'
include { GENEPREDICTION } from '../subworkflows/local/GENE_PREDICTION'
include { BINNING } from '../subworkflows/local/BINNING'
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

//SVFLOWINPUT subworkflow input
ch_contig_info = Channel.empty() //CONTIGSTAT.out.contig_info || ASSEMBLY.out.contig_info
ch_marker_profile = Channel.empty() //RAPID_TAXONOMIC_PROFILING.out.ch_profile
ch_geneset_profile = Channel.empty() //GENESET.out.ch_profile
ch_bins_list = Channel.empty() //BINNER.out.bins_list
ch_bins_folder = Channel.empty() //BINNER.out.bins_folder
ch_bins_count = Channel.empty() //BINNER.out.bins_count
ch_qs_quality_report = Channel.empty() //BINNER.out.qs_quality_report
ch_bins_info = Channel.empty() //BINNER.out.bins_info 
ch_bins_rename_map = Channel.empty() //BINNER.out.bins_rename_map
ch_gtdb_result = Channel.empty() //BINTAXONOMY.out.gtdb_summary
ch_depth_list = Channel.empty() //BINABUNDANCE.out.depth_list
ch_bins_rel_abun = Channel.empty() //BINABUNDANCE.out.totalRelativeAbun
ch_bins_count_abun = Channel.empty() //BINABUNDANCE.out.totalCountAbun
ch_bins_mean_abun = Channel.empty() //BINABUNDANCE.out.totalMeanAbun

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
            INPUT_CHECK("raw", sample_number, ch_input)
        }else{ exit 1, "Invalid input samplesheet: expects column raw_reads1 and raw_reads2 or raw_se !"}
    //Modes 2 and 3 can accept raw reads.
    }else if( params.mode in [2, 3] ){
        if(columns.contains("clean")) {
            INPUT_CHECK("clean", sample_number, ch_input)
        }else if(columns.contains("raw")){
            INPUT_CHECK("raw", sample_number, ch_input)
            raw_reads_flag = true
        }else{ exit 1, "Invalid input samplesheet: expects column raw_reads1,raw_reads2,raw_se or clean_reads1,clean_read2,clean_se !"}
    //For inputting clean reads.
    }else{
        // println("clean_reads")
        if(columns.contains("clean")) { 
            INPUT_CHECK("clean", sample_number, ch_input)
        }else{ exit 1, "Invalid input samplesheet: expects column clean_reads !"}
    }
    
    ch_raw_reads = INPUT_CHECK.out.raw_reads

    ch_clean_reads = INPUT_CHECK.out.clean_reads        
    
    CONTIGFILTER(INPUT_CHECK.out.contig)
    ch_contig = CONTIGFILTER.out.contigs

    //Submodule Execution.
    switch(params.mode){
        //QC
        case 1: {
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
            // PRODIGAL (ch_contig)

            def require_prodigal = true
            if(params.HQ_unique_bins || params.rawbin_info){
                require_prodigal = false
            }

            if(params.bin_function){
                //Check Input Parameters.
                extraParamsList = [ params.gene_info, params.cdhit_clstr, params.emapper_annotation ]
                extra_flag = checkEssentialParams(extraParamsList)
                
                //When the required configuration files (three) exist, only execute prodigal; otherwise, execute the GENESET submodule.
                if(extra_flag){
                    //No bins available then use protein for binning.
                    if(require_prodigal && !params.prodigal_output){ exit 1, "Error: the parameter --prodigal_output is required !"}
                    ch_gene_info = params2Channel(params.gene_info)
                    ch_cdhit_clstr = params2Channel(params.cdhit_clstr)
                    ch_annotation = params2Channel(params.emapper_annotation)
                    ch_vfdb_anno = params2Channel(params.vfdb_annotation)
                    ch_card_anno = params2Channel(params.card_annotation)
                }else{
                    
                    if(params.HQ_unique_bins || params.rawbin_info){
                        if(!params.prodigal_output){ exit 1, "Error: the parameter --prodigal_output is required !"}
                    }

                    require_prodigal = false
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
            }

            if(require_prodigal){
                GENEPREDICTION(sample_number, ch_contig)
                ch_prodigal_faa = GENEPREDICTION.out.prodigal_faa
                ch_report_input = ch_report_input.mix(GENEPREDICTION.out.prodigal_log)
            }

            BINNING(sample_number, ch_contig, ch_clean_reads, ch_input, ch_prodigal_faa, ch_gene_info, ch_cdhit_clstr, ch_annotation, ch_vfdb_anno, ch_card_anno)
            ch_bins_list = BINNING.out.bins_list
            ch_bins_folder = BINNING.out.bins_folder
            ch_bins_count = BINNING.out.bins_count
            ch_qs_quality_report = BINNING.out.qs_quality_report
            ch_bins_info = BINNING.out.bins_info
            ch_bins_rename_map = BINNING.out.bins_rename_map
            ch_gtdb_result = BINNING.out.gtdb_result
            ch_depth_list = BINNING.out.depth_list
            ch_bins_rel_abun = BINNING.out.bins_rel_abun
            ch_report_input = ch_report_input.mix(BINNING.out.binning_report)
            ch_bowtie2_log = ch_bowtie2_log.mix(BINNING.out.binning_bowtie2)

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
                BINNING(sample_number, ch_contig, ch_clean_reads, ch_input, ch_prodigal_faa, ch_gene_info, ch_cdhit_clstr, ch_annotation, ch_vfdb_anno, ch_card_anno)
                ch_bins_list = BINNING.out.bins_list
                ch_bins_folder = BINNING.out.bins_folder
                ch_bins_count = BINNING.out.bins_count
                ch_qs_quality_report = BINNING.out.qs_quality_report
                ch_bins_info = BINNING.out.bins_info
                ch_bins_rename_map = BINNING.out.bins_rename_map
                ch_gtdb_result = BINNING.out.gtdb_result
                ch_depth_list = BINNING.out.depth_list
                ch_bins_rel_abun = BINNING.out.bins_rel_abun
                ch_report_input = ch_report_input.mix(BINNING.out.binning_report)
                ch_bowtie2_log = ch_bowtie2_log.mix(BINNING.out.binning_bowtie2)
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
