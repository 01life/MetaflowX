
process MULTIBINNERWARNING {
    
    tag "$id"
    
    label 'process_low'

    input:
    
    tuple val(id),val(binner_number),val(tsv_number)

    output:
    path("${id}_binner_warning.log"),emit:"log"

    when:
    binner_number != tsv_number

    script:
    """
    cat <<-OUTLOG > ${id}_binner_warning.log

    ==========Start at : `date` ==========
    ### Step ${task.process} 
    ### Sample ${id}
    The total number of executions of the binning tool is: ${binner_number}, and only ${tsv_number} binning results are available. Please check for any binners that did not execute successfully. If excluding non-input file and software environment issues, it is highly likely that this binner is not performing well on this data, unable to generate high-quality bins. It is recommended that you replace the binner or exclude this sample.
    ==========End at : `date` ==========

    OUTLOG

    """
    stub:
    """
    touch ${id}_binner_warning.log
    """
}
