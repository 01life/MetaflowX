process MERGEREBINUNIMPROVEBIN {

    label 'process_single'

    input:
    path(unimprovefa)
    
    output:
    path("ReAss_ReBin_Unimprove_Bin_ID.txt"),emit:"unimprove"

    when:
    task.ext.when == null || task.ext.when

    script:
   
    """

    ls ${unimprovefa} | sed 's/^ReAss_ReBin_Unimprove_//g' |sed 's/.fa//g' > ReAss_ReBin_Unimprove_Bin_ID.txt

    """
}
