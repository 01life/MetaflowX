process PICKREBIN {

    tag "$id"

    label 'process_low'
    
    input:
    tuple val(id),path(checkm2_quality_report),path(das_bins),path(bin_QS_taxonomy)
    

    output:
    tuple val(id),path("ReAss_ReBin_*fa"),emit:"rebinfa"
    path("${id}_ReAss_ReBin_checkm2.txt"),emit:"rebinCheckm2"
    path("${id}_ReBin_checkm2.txt"),emit:"rebin_org_Checkm2"
    path("${id}_ReAss_ReBin_checkm2.txt_improve_info.txt"),emit:"improved_info"

    when:
    task.ext.when == null || task.ext.when

    script:
   
    """
    bra_pick_HQ_rebinning_bin_dastools.py --rebin_checkm ${checkm2_quality_report}  --das_tools_result ${das_bins} --org_checkm  ${bin_QS_taxonomy} --bin_id ${id} --mergeQS ${id}_ReAss_ReBin_checkm2.txt  --pick_checkm ${id}_ReBin_checkm2.txt
    """

    stub:
    """
    touch ReAss_ReBin_0.fa
    touch ${id}_ReAss_ReBin_checkm2.txt
    touch ${id}_ReBin_checkm2.txt
    touch ${id}_ReAss_ReBin_checkm2.txt_improve_info.txt
    """

}
