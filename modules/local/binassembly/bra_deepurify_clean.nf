process DEEPURIFYCLEAN {

    tag "$id"

    label 'process_high'

    input:
    tuple val(id),path(refine_bin)
    path(checkm2_db)

    output:
    tuple val(id),path("${id}_deepurify_bin"),emit:"deepurifyDir"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.deepurify_clean_options ?: ""

    """
    export CHECKM2DB="${checkm2_db}"

    if [ -f "${refine_bin}" ]; then
        mkdir -p rowBin
        mv "${refine_bin}" rowBin/
    elif [ -d "${refine_bin}" ]; then
        ln -s "${refine_bin}" rowBin
    else
        echo "Error: ${refine_bin} is neither a file nor a directory."
        exit 1
    fi

    #version 2.3.3
    
    deepurify clean  \\
        --num_process ${task.cpus} \\
        --db_folder_path ${params.deepurify_db} \\
        -i rowBin  \\
        -o ${id}_deepurify_bin/ \\
        ${options}
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Deepurify: \$(echo \$(deepurify --version 2>&1) | sed 's/^.*deepurify //; s/ .*\$//')
    END_VERSIONS


    """
    stub:
    """
    mkdir ${id}_deepurify_bin
    """


}

process DEEPURIFYCLEANRENAME {

    tag "$id"

    label 'process_low'

    input:
    tuple val(id),path(deepurify_raw)

    output:
    tuple val(id),path("${id}_Deepurify_rename_bin"),emit:"renameDeepurify"

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    bra_deepurify_rename_bin.py --deepurify_dir ${deepurify_raw} --rename_dir ${id}_Deepurify_rename_bin  --metaInfo  ${deepurify_raw}/MetaInfo.tsv --binid ${id}
    """

    stub:
    """
    mkdir ${id}_Deepurify_rename_bin
    """
}
