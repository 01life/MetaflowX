process PICKREREFINE {


    label 'process_low'

    input:
    path(renameDeepurify)
    path(rebining_checkm2)
    

    output:
    path("ReAss_ReBin_ReRefine_*.txt") ,emit:"rebinfa",optional: true
    path("ReAss_ReBin_ReRefine*.fa"),emit:"HQ_Rerefinefa",optional: true
    path("ReAss_ReBin_ReRefine_checkm2.txt"),emit:"refine_checkm2"
    path("ReRefine_improve_info.txt"),emit:"improved_info"

    when:
    task.ext.when == null || task.ext.when

    script:
   
    """
    mkdir merger_floder
    cp -rf ${renameDeepurify}/ ./merger_floder
    mkdir tmp
    cp -rf ./merger_floder/*/*fa ./tmp/
    cat ./merger_floder/*/deepurify_rename.QS.txt  > ReAss_ReBin_ReRefine_checkm2.txt

    bra_pick_HQ_re-refine_bin_deepurify.py --mergeQS ${rebining_checkm2} --deepurify_info ReAss_ReBin_ReRefine_checkm2.txt --deepurify_fa tmp --failed ReAss_ReBin_ReRefine_Unimprove_Bin_ID.txt --improve ReRefine_improve_info.txt
    """

    stub:
    """
    touch ReAss_ReBin_ReRefine_checkm2.txt
    touch ReRefine_improve_info.txt
    touch ReAss_ReBin_ReRefine_Unimprove_Bin_ID.txt
    touch ReAss_ReBin_ReRefine_1.txt
    touch ReAss_ReBin_ReRefine1.fa
    """
}
