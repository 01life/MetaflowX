
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { checkEssentialParams } from '../modules/local/common/utils'
include { checkFaFiles } from '../modules/local/common/utils'

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetassembly.initialise(params, log)


if(params.bin_genome_floder) { checkFaFiles(params.bin_genome_floder) }

if (params.outdir) { ch_output = new File(params.outdir).getAbsolutePath()  } else { exit 1, "Output directory not specified !" }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { EASYREASSEMBLY } from '../subworkflows/local/EASY_REASSEMBLY'
include { INPUT_CHECK } from '../subworkflows/local/INPUT_CHECK'

include { BINMAP } from '../subworkflows/local/BIN_MAP'
// include { BINREFINE } from '../subworkflows/local/BIN_REFINE'
include { BINREASSEMBLY} from '../subworkflows/local/BIN_REASSEMBLY'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow CUSTOMREASSEMBLY {

    bra_essential_db = [params.input, params.bin_QS_taxonomy, params.bins_count_abun, params.bins_mean_abun, params.bin_genome_floder]

    if(!checkEssentialParams(bra_essential_db)) { exit 1, "The required parameters to execute the BinReAssembly module  are:\n --input\n --bin_QS_taxonomy\n --bins_count_abun \n --bins_mean_abun\n --bin_genome_floder" }


    //Input Check and Data Parsing
    if (params.input) { ch_input = file(params.input) } else { exit 1, "Input samplesheet not specified !" }

    
    //Determine if parsing is for raw_reads/clean_reads based on column names.
    def columns = ch_input.text.tokenize('\n').first()

    if(columns.contains("clean")) { 
        INPUT_CHECK("clean", ch_input)
    }else{ exit 1, "Invalid input samplesheet: expects column clean_reads !"}

    if (params.outdir) { ch_output = new File(params.outdir).getAbsolutePath()  } else { exit 1, "Output directory not specified !" }


    ch_clean_reads = INPUT_CHECK.out.clean_reads 

    ch_bin_genome = Channel.fromPath("${params.bin_genome_floder}/*.fa")
    ch_bin_QS_taxonomy = Channel.fromPath(params.bin_QS_taxonomy) 
    ch_bins_count_abun = Channel.fromPath(params.bins_count_abun)
    ch_bins_mean_abun = Channel.fromPath(params.bins_mean_abun)

    BINMAP(ch_clean_reads, ch_bin_QS_taxonomy, ch_bins_count_abun, ch_bins_mean_abun, ch_bin_genome)
    ch_pick_info = BINMAP.out.ch_bin_sample_ref_reads_binfa
    ch_binsample = BINMAP.out.binsample

    // bin refine
    // if(params.bin_refine){
    //     BINREFINE(ch_output, ch_clean_reads, ch_contig, ch_filter_bin, ch_filter_bin_info, ch_filter_bin_mash_fq)
    //     ch_filter_bin_mash_fq = BINREFINE.out.after_refine_bin_mash_fq
    // }

    //bin reassembly

    BINREASSEMBLY(ch_pick_info, ch_binsample, ch_clean_reads, ch_bin_QS_taxonomy)
    
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
