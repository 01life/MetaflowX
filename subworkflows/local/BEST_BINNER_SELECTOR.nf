include { COMBINEBINNER } from '../../modules/local/bbs/combineBinner2DASTools'
include { MULTIDASTOOL } from '../../modules/local/bbs/multi_DASTools'
include { PIPELINEWARNING } from '../../modules/local/common/pipeline_warning'
include { CATEVALBYID } from '../../modules/local/bbs/cat_eval'
include { SUMMARYRESULT } from '../../modules/local/bbs/summary_DASTools'
include { ZIPDASTOOLRES } from '../../modules/local/bbs/zip_dastool_res'


workflow BEST_BINNER_SELECTOR {
    take:
    contig2bin
    ch_protein
    ch_contig

    main:

    COMBINEBINNER(contig2bin)

    ch_eval_tsv = COMBINEBINNER.out.binner_combination
        .flatten()
        .map { it -> [it.baseName, it] }

    //path(dastoolinput),path(contig),path(protein)
    ch_dastool = ch_eval_tsv.join(ch_protein).join(ch_contig)

    MULTIDASTOOL(ch_dastool)
    PIPELINEWARNING("BBS_DASTool", MULTIDASTOOL.out.das_bins_error.collect())

    CATEVALBYID(MULTIDASTOOL.out.eval)
    SUMMARYRESULT(CATEVALBYID.out.eval.map{ it -> it[1] }.collectFile(name: 'all_eval_noheader.xls'))

    ch_eval_tsv = MULTIDASTOOL.out.eval.join(MULTIDASTOOL.out.tsv)
    ZIPDASTOOLRES(ch_eval_tsv)

}

