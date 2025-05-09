process MERGECHECKLOG {

    label 'process_low'

    input:
    val(sample_number)
    val(finish_number)
    path(check_log)

    output:
    path("${params.pipeline_prefix}_INPUT_CHECK_log_info.txt"),emit:"info"

    when:
    sample_number == finish_number

    script:
    """

    cat ${check_log} > ${params.pipeline_prefix}_INPUT_CHECK_log_info.txt

    """

}