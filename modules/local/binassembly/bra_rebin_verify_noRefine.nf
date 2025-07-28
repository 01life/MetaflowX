
process REBINVERIFY {

    label 'process_low'

    input:
    path(pickbininfo)
    path(orignal)
    path(reass)
    path(rebin)
    path(refine)
    path(reass_improved_info)
    path(rebin_improved_info)
    path(refine_improved_info)

    output:
    path("bin_res_report"),emit:"bin_reassembly_stat_info"
    path("Bins_Reassembly_Optimization_evaluation_info.xls")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${reass_improved_info} ${rebin_improved_info} ${refine_improved_info} > merged_HQ_improve.txt

    bra_reassembly_bins_quality_evaluator.py  --pickbin  ${pickbininfo} --original_QS  ${orignal} --ReAss_QS   ${reass} --ReBin_QS   ${rebin} --ReRefine_QS  ${refine} --outfile Bins_Reassembly_Optimization_evaluation_info.xls --platout Bins_Reassembly_Optimization_evaluation --improve  merged_HQ_improve.txt


    mkdir bin_res_report 
    cp Bins_Reassembly_Optimization_evaluation_final_improvement.xls bin_res_report/bin_res_stat.xls
    cp Bins_Reassembly_Optimization_evaluation.point.txt bin_res_report/bin_res.point.xls
    cp Bins_Reassembly_Optimization_evaluation.sanket.txt bin_res_report/bin_res.sanket.xls

    """

    stub:
    """
    touch Bins_Reassembly_Optimization_evaluation_info.xls
    mkdir bin_res_report
    """

}
