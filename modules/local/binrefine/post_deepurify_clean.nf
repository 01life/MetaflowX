process POSTDEEPURIFYCLEAN {

    label 'process_low'

    publishDir(
        path: "${params.outdir}/08.BinOptimization/081.BinRefine/Deepurify/",
        pattern: "Deepurify_Result",
        mode: "copy",
        failOnError: true
    )

    input:
    path(refine_bin)
    path(filter_bin_mash_fq)
    path(deepurify_bin)


    output:
    path("Deepurify_rename_bin/*.fa"),emit:"deepurify_bin"
    path("${params.pipeline_prefix}_deepurify_bin_mash_fq.txt"),emit:"deepurify_bin_mash_fq"
    path("Deepurify_Result"),emit:"deepurify_raw"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    deepurify_rename_bin.py ${refine_bin} ${deepurify_bin} Deepurify_rename_bin ${deepurify_bin}/MetaInfo.tsv 


    replace_deepurify_bin_2_mash_fq_map.py ./Deepurify_rename_bin/ ${filter_bin_mash_fq} ${params.pipeline_prefix}_deepurify_bin_mash_fq.txt

    cp -rf ${deepurify_bin}  Deepurify_Result

    mkdir -p Deepurify_Result/rename_bins/

    cp -rf Deepurify_rename_bin/*.fa Deepurify_Result/rename_bins/
    cp -rf Deepurify_rename_bin/deepurify_rename*txt Deepurify_Result/rename_bins/


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(echo \$(python --version 2>&1) | sed 's/^.*python //; s/ .*\$//')
    END_VERSIONS


    """

}
