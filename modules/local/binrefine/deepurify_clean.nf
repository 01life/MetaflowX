process DEEPURIFYCLEAN {

    tag "$id"

    label 'process_high'

    input:
    tuple val(id),path(refine_bin)
    path(checkm2_db)

    output:
    tuple val(id),path("${id}_deepurify_bin"),emit:"res"
    tuple val(id),path("${id}_Deepurify_rename_bin"),emit:"renameDeepurify"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.deepurify_clean_options ?: ""
    """
    export CHECKM2DB="${checkm2_db}"

    if [ -f "${refine_bin}" ]; then
        mkdir -p tmp
        mv "${refine_bin}" tmp/
    elif [ -d "${refine_bin}" ]; then
        ln -s "${refine_bin}" tmp
    else
        echo "Error: ${refine_bin} is neither a file nor a directory."
        exit 1
    fi

    #version 2.3.3
    
    deepurify clean  \\
        --num_process ${task.cpus} \\
        --db_folder_path ${params.deepurify_db} \\
        -i tmp  \\
        -o ${id}_deepurify_bin/ \\
        ${options}

    deepurify_rename_bin.py tmp ${id}_deepurify_bin ${id}_Deepurify_rename_bin ${id}_deepurify_bin/MetaInfo.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Deepurify: \$(echo \$(deepurify --version 2>&1) | sed 's/^.*deepurify //; s/ .*\$//')
    END_VERSIONS


    """

}
