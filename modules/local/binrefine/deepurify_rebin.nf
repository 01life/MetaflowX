process DEEPURIFYREBIN {

    label 'process_high'

    input:
    path(contigs)
    path(sorted_bam)
    path(checkm2_db)

    output:
    path("deepurify_bin"),emit:"res"
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.deepurify_rebin_options ?: ""
    """
    export CHECKM2DB="${checkm2_db}"

    #version 2.3.3
    deepurify re-bin  \\
        --num_process ${task.cpus} \\
        --db_folder_path ${params.deepurify_db} \\
        -c ${contigs}  \\
        -s ${sorted_bam} \\
        -o deepurify_bin/ \\
        ${options}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Deepurify: \$(echo \$(deepurify --version 2>&1) | sed 's/^.*deepurify //; s/ .*\$//')
    END_VERSIONS


    """

}
