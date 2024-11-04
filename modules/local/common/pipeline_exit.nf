
process PIPELINEEXIT {

    label 'process_low','error_finish'

    input:
    path(error_log)

    when:
    task.ext.when == null || task.ext.when

    script:
    def logList = error_log instanceof List ? error_log.join(",") : "$error_log"
    """
    
    echo "The pipeline execution encountered an error, please refer to the log file: ${logList}."

    exit 1

    """

}
