
process FILTERREFINE {

    label 'process_low'
    
    publishDir(
        path: "${params.outdir}/08.BinOptimization/",
        pattern: "*Refine_bin*info.txt",
        mode: "copy",
        failOnError: true
    )

    input:
    path(priginal_bin_QS_taxonomy)
    path(refine_bin_QS_taxonomy)
    path(refine_bin_list)
    path(reassembly_map)


    output:

    path("after_refine_bin_mash_fq.txt"),emit:"after_refine_bin_mash_fq"
    path("Deepurify_COBRA_Refine_bin_choose_info.txt"),emit:"refine_map"
    path("After_Refine_bin_CheckM2_info.txt"),emit:"refine_checkm2"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    verify_refine_result.py ${priginal_bin_QS_taxonomy} ${refine_bin_QS_taxonomy} ${refine_bin_list} ${reassembly_map} after_refine_bin_mash_fq.txt

    cp ${refine_bin_QS_taxonomy} After_Refine_bin_CheckM2_info.txt
    cp ${refine_bin_list} Deepurify_COBRA_Refine_bin_choose_info.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS


    """

}
