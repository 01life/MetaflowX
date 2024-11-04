
process PIPELINEERROR {

    label 'process_low'

    input:
    val(step_name)
    val(sample_number)
    val(finish_number)

    output:
    path("${params.pipeline_prefix}_${step_name}_error_log*.txt"),emit:"log"

    when:
    sample_number != finish_number

    script:
    """

    timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')  

    cat <<-OUTLOG > ${params.pipeline_prefix}_${step_name}_error_log_\$timestamp.txt

    ==========Start at : `date` ==========
    ### Step ${task.process}
    There are ${sample_number} samples in total, with ${finish_number} samples successfully completing the ${step_name} step.
    ==========End at : `date` ==========

    OUTLOG

    """

}
