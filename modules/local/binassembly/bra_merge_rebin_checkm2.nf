libgcc=process MERGEREBINCHECKM2 {

    label 'process_single'
    
    input:
    path(rebinCheckm2)
    path(improve)

    output:

    path("ReAss_ReBin_checkm2.txt"),emit:"rebining_checkm2"
    path("ReBin_improve_info.txt"),emit:"improved_info"

    when:
    task.ext.when == null || task.ext.when

    script:
   
    """

    # Get the first filename in the rebinCheckm2 directory and write its first line to c1.
    head -n 1 "\$(ls ${rebinCheckm2} | head -n 1)" > c1

    # Get the contents of all files in the rebinCheckm2 directory (starting from the second line) and write to c2.
    for file in ${rebinCheckm2}; do
        tail -n +2 "\$file"
    done > c2

    # Merge c1 and c2 into ReAss_ReBin_checkm2.txt.
    cat c1 c2 > ReAss_ReBin_checkm2.txt

    # Delete temporary files c1 and c2.
    rm c1 c2


    cat ${improve} > ReBin_improve_info.txt


    """
}
