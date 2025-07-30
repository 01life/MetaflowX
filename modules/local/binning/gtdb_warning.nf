process GTDBWARNING {
    
    label 'process_low'

    input:
    val(summary_linecount)

    output:
    path("${params.pipeline_prefix}_GTDBTk_warning_log*.txt"), emit: "log"

    when:
    summary_linecount < 2

    script:
    """
    timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')  
    
    cat <<-OUTLOG > ${params.pipeline_prefix}_GTDBTk_warning_log_\$timestamp.txt
    
    ==========Start at : `date` ==========
    ### Step ${task.process}
    Do not get gtdbtk.bac120.summary.tsv and gtdbtk.ar53.summary.tsv, please cheak your result!
    ==========End at : `date` ==========

    OUTLOG

    """
    
    stub:
    """
    timestamp=\$(date '+%Y-%m-%d_%H-%M-%S')  
    touch ${params.pipeline_prefix}_GTDBTk_warning_log_\$timestamp.txt
    """
}
