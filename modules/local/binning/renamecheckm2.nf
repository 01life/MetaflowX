
process RENAMECHECKM2 {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(quality_report)
    val(step_name)

    output:
    tuple val(id), path("${id}_${step_name}_quality_report.tsv"), emit: "new_quality_report"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    cp ${quality_report} ${id}_${step_name}_quality_report.tsv
    
    """
    stub:
    """
    touch ${id}_${step_name}_quality_report.tsv
    """
}
