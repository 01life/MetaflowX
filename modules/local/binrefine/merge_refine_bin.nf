
process MERGEREFINEBIN {

    label 'process_single'

    publishDir(
        path: "${params.outdir}/08.BinOptimization/081.BinRefine/Last_Refined_Bin/",
        pattern: "Refine_bin",
        mode: "copy",
        failOnError: true
    )


    input:
    path(deepurify_bin)
    path(cobra_bin)
    
    output:
    tuple val("after_refine_bin"),path("Refine_bin"),emit:"refine_bin"
    path("refinebin_path.txt"),emit:"refine_bin_list"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    mkdir deepurify 
    cp ${deepurify_bin} ./deepurify

    mkdir cobra 
    cp ${cobra_bin} ./cobra

    verify_deepurify_cobra_bin.py  ./deepurify  ./cobra ./Refine_bin refinebin_path.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS


    """

}
