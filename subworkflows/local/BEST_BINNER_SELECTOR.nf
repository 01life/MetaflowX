include { COMBINEBINNER } from '../../modules/local/bbs/combineBinner2DASTools'
include { MULTIDASTOOL } from '../../modules/local/bbs/multi_DASTools'
include { PIPELINEWARNING } from '../../modules/local/common/pipeline_warning'
include { SUMMARYRESULT } from '../../modules/local/bbs/summary_DASTools'
include { ZIPDASTOOLRES } from '../../modules/local/bbs/zip_dastool_res'


workflow BEST_BINNER_SELECTOR {
    take:
    contig2bin
    protein_list
    contig_list

    main:

    COMBINEBINNER(contig2bin, protein_list, contig_list)

    ch_dastool = COMBINEBINNER.out.binner_combination.splitText()
        .map { it -> 
            def split = it.trim().split("\t")
            def id = split[0]
            def contig = split[1]
            def protein = split[2]
            def label = split[3]
            // collect(): Iterate through List.
            def contig2bin = split[-1].trim().split(",").collect { file(it) }
            return [ id, contig, protein, label, contig2bin ]
        }
    MULTIDASTOOL(ch_dastool)
    PIPELINEWARNING("BBS_DASTool", MULTIDASTOOL.out.das_bins_error.collect())

    SUMMARYRESULT(MULTIDASTOOL.out.eval.collect())

    ZIPDASTOOLRES(MULTIDASTOOL.out.eval.collect(), MULTIDASTOOL.out.tsv.collect(), contig_list)

}

