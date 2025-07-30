process NOBINSWARNING {
    
    label 'process_low'

    input:
    val(filtered_pass)

    output:
    path("${params.pipeline_prefix}_Binning_warning_log*.txt"), emit: "log"

    when:
    filtered_pass == 0

    script:
    """
    timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')  
    
    cat <<-OUTLOG > ${params.pipeline_prefix}_Binning_warning_log_\$timestamp.txt
    
    ==========Start at : `date` ==========
    ### Step ${task.process}
    Oppo~ Do not get any bins in all samples. Please check your input and data. May be your data do not suit by asseembly. SO BAD~~~ T -_- T
    ==========End at : `date` ==========

    OUTLOG

    """

    stub:
    """
    timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')
    cat <<-OUTLOG > ${params.pipeline_prefix}_Binning_warning_log_\$timestamp.txt
    ==========Start at : `date` ==========
    ### Step ${task.process}
    Oppo~ Do not get any bins in all samples. Please check your input and data. May be your data do not suit by asseembly. SO BAD~~~ T -_- T
    ==========End at : `date` ==========
    OUTLOG
    """
    
}
