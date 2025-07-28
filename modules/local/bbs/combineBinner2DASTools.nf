process COMBINEBINNER {
    
    label 'process_low'

    input:
    path(contig2bin)

    output:
    path("AllSample/*"),emit:"binner_combination"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    create_multiDASTools_input.py -t ${contig2bin} -o AllSample

    """
    stub:
    """
    mkdir -p AllSample/Sample1
    mkdir -p AllSample/Sample2
    mkdir -p AllSample/Sample3
    """
    
}