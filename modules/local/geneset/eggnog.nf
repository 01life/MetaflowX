
process EGGNOG {
    
    tag "$basename"

    label 'process_medium'

    input:
    tuple val(basename),path(pep)
    path(eggnog_diamond_db)
    path(eggnog_mapper_db)

    output:
    path("${basename}.emapper.annotations"),emit:"annotation"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.eggnog_options ?: ""
    
    """

    emapper.py -i ${pep} -o ${basename} --dmnd_db ${eggnog_diamond_db} --data_dir ${eggnog_mapper_db} --cpu ${task.cpus} ${options}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(echo \$(emapper.py --version) | grep -o "emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+" | sed "s/emapper-//")
    END_VERSIONS

    """

    stub:
    """
    mkdir -p ${basename}.emapper.annotations
    touch ${basename}.emapper.annotations/annotations.tsv
    touch ${basename}.emapper.annotations/names.tsv
    """
    
}
