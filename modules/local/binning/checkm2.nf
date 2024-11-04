
process CHECKM2 {
    
    tag "$id"
    
    label 'process_single'

    input:
    tuple val(id),path(bins)
    path(checkm2_db)

    output:
    tuple val(id),path("${id}/quality_report.tsv"), emit: "quality_report", optional: true
    path("${id}_checkm2_error.txt"), emit: "checkm2_error", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.checkm2_options ?: ""

    """

    checkm2 predict --threads ${task.cpus} --input ${bins} --output-directory ${id} --database_path ${checkm2_db} ${options} || {
    
    cat <<-OUTLOG > ${id}_checkm2_error.txt

    ==========Start at : `date` ==========
    ### Step ${task.process}
    CheckM2 task for sample ${id} failed. Please check if bins are present in folder ${bins}.
    ==========End at : `date` ==========

    OUTLOG
    
    }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
        
    """

}
