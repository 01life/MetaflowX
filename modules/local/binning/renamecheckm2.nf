
process RENAMECHECKM2 {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(quality_report)

    output:
    tuple val(id), path("${id}_quality_report.tsv"), emit: "new_quality_report"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    cp ${quality_report} ${id}_quality_report.tsv
    
    """
}
