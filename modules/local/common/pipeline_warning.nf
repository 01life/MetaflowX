process PIPELINEWARNING {
    
    label 'process_low'

    input:
    val(step_name)
    path(warning_log)

    output:
    path("${params.pipeline_prefix}_${step_name}_warning_log*.txt"), emit: "log"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')  
    cat ${warning_log} > ${params.pipeline_prefix}_${step_name}_warning_log_\$timestamp.txt

    """

}
