//
// Assembly
//


include { MEGAHIT } from '../../modules/local/assembly/megahit'
include { SWITCH2MEGAHIT } from '../../modules/local/assembly/switch2megahit'
include { METASPADES } from '../../modules/local/assembly/metaspades'
include { METASPADESM64 } from '../../modules/local/assembly/metaspades_m64'
include { METASPADESM128 } from '../../modules/local/assembly/metaspades_m128'
include { METAQUAST } from '../../modules/local/assembly/metaquast'
include { MERGEQUAST } from '../../modules/local/assembly/merge_metaquast'
include { CONTIGSTAT } from '../../modules/local/assembly/contig_stat'
include { PIPELINEWARNING } from '../../modules/local/common/pipeline_warning'
include { PIPELINEERROR } from '../../modules/local/common/pipeline_error'
include { PIPELINEEXIT } from '../../modules/local/common/pipeline_exit'


workflow ASSEMBLY {
    take:
    sample_number
    clean_reads      // channel: [ val(id), [ reads1, reads2 ] ]

    main:

    if (!(params.assembly_tool in ["metaspades", "megahit"])) { exit 1, "The parameter assembly_tool is invalid, supported values are:\n * metaspades \n * megahit" }

    ch_warning_log = Channel.empty()

    if (params.assembly_tool == "metaspades") {
        
        // METASPADES (clean_reads)
        // spades_contigs = METASPADES.out.contigs

        //metaspades q_32_64
        METASPADESM64(clean_reads)
        spades_m64_contigs = METASPADESM64.out.contigs

        //metaspades q_32_64 failed => q_32_128
        ch_spades_m64_failed = clean_reads.join(spades_m64_contigs, remainder: true).filter{ it[2] == null }.map{ it -> [it[0], it[1]]}
        METASPADESM128(ch_spades_m64_failed)
        spades_m128_contigs = METASPADESM128.out.contigs

        // get samples which MetaSPAdes completed
        spades_contigs = spades_m64_contigs.concat(spades_m128_contigs).filter{ id, contigs -> contigs.size() > 0 }

        // get samples which MetaSPAdes failed and switch to MegaHit for assembly
        ch_spades_failed = clean_reads.join(spades_contigs, remainder: true).filter{ it[2] == null }.map{ it -> [it[0], it[1]]}
        MEGAHIT(ch_spades_failed)
        megahit_contigs = MEGAHIT.out.contigs

        SWITCH2MEGAHIT(megahit_contigs)
        // merge and publish warning log
        PIPELINEWARNING("Assembly", SWITCH2MEGAHIT.out.log.collect())
        ch_warning_log = PIPELINEWARNING.out.log

        megahit_contigs = megahit_contigs.filter{ id, contigs -> contigs.size() > 0 }
        contigs = spades_contigs.concat(megahit_contigs)
        
    }
    
    if (params.assembly_tool == "megahit") {
        MEGAHIT (clean_reads)
        contigs = MEGAHIT.out.contigs
    }

    all_contig = contigs.map{ _,file -> file}.collect()
    finish_number = all_contig.flatten().filter { file -> !file.isEmpty() }.count()
    
    // assembly completed successfully
    // run the assembly validation step
    METAQUAST (contigs)
    MERGEQUAST(METAQUAST.out.report.collect())

    CONTIGSTAT(sample_number, finish_number, all_contig)
    contig_info = CONTIGSTAT.out.contig_info
    contig_report = CONTIGSTAT.out.contig_report

    // generate an error log and terminate the pipeline if Assembly error occurs
    PIPELINEERROR("Assembly", sample_number, finish_number)
    PIPELINEEXIT(PIPELINEERROR.out.log)

    assembly_report = contig_report.mix(ch_warning_log, PIPELINEERROR.out.log).collect()

    emit:
    contigs           // channel: [ val(id), [ contigs ] ]
    assembly_report
    contig_info
    
}