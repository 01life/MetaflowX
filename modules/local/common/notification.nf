
process NOTIFICATION {

    label 'process_low'

    input:
    val(step_name)
    val(content)
    path(notification)
    
    output:
    path("${params.pipeline_prefix}_${step_name}_notification_log.txt"),emit:"log"

    when:
    task.ext.when == null || task.ext.when

    script:
    def webhookurl = params.webhookurl ?: ""
    def username = System.getProperty('user.name')
    """

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

    cp ${notification} ${params.pipeline_prefix}_${step_name}_notification_log.txt

    """
    
}
