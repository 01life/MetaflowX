//
// geneset construction
//

include { GENEFILTERSTAT } from '../../modules/local/geneset/gene_filter_stat'
include { MERGEGENEINFO } from '../../modules/local/geneset/merge_gene_info'
include { SINGLECDHIT } from '../../modules/local/geneset/single_cdhit'
include { MULTICDHIT } from '../../modules/local/geneset/multi_cdhit'
include { CDHITDIV } from '../../modules/local/geneset/multi_cdhit_1div'
include { ESTDIV } from '../../modules/local/geneset/multi_cdhit_2est_div'
include { CDHITEST2D } from '../../modules/local/geneset/multi_cdhit_3est2d'
include { CDHITCLSTR } from '../../modules/local/geneset/multi_cdhit_4create_mergeclstr'
include { MEGRCLSTR } from '../../modules/local/geneset/multi_cdhit_5mergeclstr'
include { PASTERESULT } from '../../modules/local/geneset/multi_cdhit_6paste_result'
// include { CDHIT } from '../../modules/local/geneset/rename_cdhit'

workflow GENESETCONSTRUCTION {
    take:
    sample_number
    prodigal_cds   
    prodigal_faa   
    // prodigal_noid_cds   //PRODIGAL.out.noid_cds.collect()
    // prodigal_noid_faa   //PRODIGAL.out.noid_faa.collect()

    main:

    prodigal_output = prodigal_cds.join(prodigal_faa)
    GENEFILTERSTAT(prodigal_output)

    all_stat_info = GENEFILTERSTAT.out.gene_stat.collectFile(name:"geneset_gene_stat.txt", keepHeader: true)

    //CD-hit support multi tasks
    prodigal_noid_faa = prodigal_faa.map{ it -> it[1] }.collect()
    finish_number_prodigal = prodigal_noid_faa.flatten().filter { file -> !file.isEmpty() }.count()

    MERGEGENEINFO(sample_number, finish_number_prodigal, all_stat_info)

    // single cd-hit task
    SINGLECDHIT(GENEFILTERSTAT.out.cds.collect(), GENEFILTERSTAT.out.pep.collect(), MERGEGENEINFO.out.single_task)

    // multi cd-hit task
    MULTICDHIT(GENEFILTERSTAT.out.cds.collect(), GENEFILTERSTAT.out.pep.collect(), MERGEGENEINFO.out.multi_task)

    // CDHITDIV(prodigal_noid_cds, GENEFILTER.out.multi_task)
    // CDHITDIV(GENEFILTER.out.allcds, GENEFILTER.out.multi_task)
    // ESTDIV(CDHITDIV.out.div.flatten())
    // ch_div = CDHITDIV.out.div.flatten().map { file ->
    //         def fileName = file.name.replaceFirst(/^all.cds.fa.div-/, '')
    //             [fileName]
    //         }
    // ch_est_input = ch_div.combine(CDHITDIV.out.div_tmp).combine(ESTDIV.out.div0.collect().toList())

    // CDHITEST2D(ch_est_input)
    // CDHITCLSTR(CDHITEST2D.out.clstr.collect(), CDHITEST2D.out.div_o.collect(), GENEFILTER.out.multi_task)
    // ch_merge_input = CDHITCLSTR.out.clstr_order
    //     .splitText()
    //     .map { it -> 
    //         def in_list = it.trim().split("\t")[1]
    //         def clstr_num = it.trim().split("\t")[0]
    //         return [clstr_num,in_list]
    //     }
    // MEGRCLSTR(ch_merge_input,CDHITEST2D.out.clstr.collect())
    // PASTERESULT(MEGRCLSTR.out.subclstr.collect(), CDHITCLSTR.out.unique_fa, GENEFILTER.out.allpep)

    // ch_cdhit_pep = SINGLECDHIT.out.pep.mix(PASTERESULT.out.pep).collect()
    // ch_cdhit_cds = SINGLECDHIT.out.cds.mix(PASTERESULT.out.cds).collect()
    // ch_cdhit_gene_length = SINGLECDHIT.out.gene_length.mix(PASTERESULT.out.gene_length).collect()
    // ch_cdhit_gene_info = SINGLECDHIT.out.gene_info.mix(PASTERESULT.out.gene_info).collect()
    // ch_cdhit_clstr = SINGLECDHIT.out.clstr.mix(PASTERESULT.out.clstr).collect()
    // ch_cdhit_split = SINGLECDHIT.out.split.mix(PASTERESULT.out.split).collect()
    // ch_cdhit_geneset_gene_report = SINGLECDHIT.out.geneset_gene_report.mix(PASTERESULT.out.geneset_gene_report).collect()

    // CDHIT(ch_cdhit_pep, ch_cdhit_cds, ch_cdhit_gene_length, ch_cdhit_gene_info, ch_cdhit_clstr, ch_cdhit_split, ch_cdhit_geneset_gene_report)
    // clstr = CDHIT.out.clstr             // channel: [ path(clstr) ]
    // split_fa = CDHIT.out.split.collect()  
    // cds = CDHIT.out.cds  
    // pep = CDHIT.out.pep
    // gene_length = CDHIT.out.gene_length
    // gene_info = CDHIT.out.gene_info     // channel: [ path(gene_info) ]
    // gene_report = CDHIT.out.geneset_gene_report

    clstr = SINGLECDHIT.out.clstr.mix(MULTICDHIT.out.clstr).collect()
    split_fa = SINGLECDHIT.out.split.mix(MULTICDHIT.out.split).collect()
    cds = SINGLECDHIT.out.cds.mix(MULTICDHIT.out.cds).collect()
    pep = SINGLECDHIT.out.pep.mix(MULTICDHIT.out.pep).collect()
    gene_info = SINGLECDHIT.out.gene_info.mix(MULTICDHIT.out.gene_info).collect()
    gene_length = SINGLECDHIT.out.gene_length.mix(MULTICDHIT.out.gene_length).collect()
    gene_report = MERGEGENEINFO.out.report.mix(SINGLECDHIT.out.geneset_gene_report, MULTICDHIT.out.geneset_gene_report).collect()

    emit:
    clstr
    split_fa
    cds
    pep
    gene_info 
    gene_length
    gene_report

}

