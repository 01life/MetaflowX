
process BIGMAPWARNING {
    
    label 'process_low'
    
    input:
    val(sample_number)
    val(finish_number)

    output:
    path("BiG-MAP_warning.log"),emit:"log"

    when:
    sample_number != finish_number

    script:
    """
    cat <<-OUTLOG > BiG-MAP_warning.log
    
    ==========Start at : `date` ==========
    ### Step ${task.process}
    The total number of samples is : ${sample_number}, and the number of samples successfully running the BiG-MAP process is ${finish_number}.
    ==========End at : `date` ==========

    OUTLOG

    """
    
}
