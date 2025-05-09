
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
    def webhookurl = params.webhookurl ?: ""
    def content = "❌ Oops! ${step_name} hit a little snag, and ${params.pipeline_prefix} is taking a break. Hurry up and check it out before it throws another tantrum! 🫨"
    def username = System.getProperty('user.name')
    """

    timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')  

    cat <<-OUTLOG > ${params.pipeline_prefix}_${step_name}_error_log_\$timestamp.txt

    ==========Start at : `date` ==========
    ### Step ${task.process}
    There are ${sample_number} samples in total, with ${finish_number} samples successfully completing the ${step_name} step.
    ==========End at : `date` ==========

    OUTLOG

    if [ -n ${webhookurl} ]; then
        curl ${webhookurl} \\
            -H 'Content-Type: application/json' \\
            -d '
            {
                "msgtype": "text",
                "text": {
                    "content": "${content}",
                    "mentioned_list": ["${username}"]
                }
            }'
    else
        echo "Friendly heads-up: If you're considering the robot alarm, make sure to add the --webhookurl parameter. If you're not interested, just ignore this message!"
    fi

    """

}
