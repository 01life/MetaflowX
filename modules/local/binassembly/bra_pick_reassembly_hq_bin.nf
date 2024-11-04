
process OUTHQBIN {

    tag "$step"

    label 'process_low'
    
    input:
    val(step)
    path(org_bin)
    path(checkm2report)

    output:
    path("01.ReAss_HQ/*fa"),emit:"hqbin",optional:true
    path("LQ_RAB/*fa"),emit:"lqbin"
    path("lq_bin_id.txt"),emit:"lqbinid"
    path("LQ_RAB"),emit:"lqbinFloder"
    path("ReAss_*.txt"),emit:"reass"
    path("ReAss_checkm2.txt"),emit:"reass_checkm2"
    path("ReAss_improved_info.txt"),emit:"improved_info"


    when:
    task.ext.when == null || task.ext.when

    script:
    def reassembly_hq_options = params.reassembly_HQ_options ?: ""
   
    """
    mkdir raw_bin
    cp ${org_bin} ./raw_bin

    bra_pick_HQ_reassembly_bin_checkm2.py --input_file ${checkm2report} --fa_dir raw_bin --hq_dir 01.ReAss_HQ --lq_bin_id_file lq_bin_id.txt ${reassembly_hq_options} 

    cp ${checkm2report} ReAss_checkm2.txt

    """
}
