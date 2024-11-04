
process MEGRCLSTR {

    label 'process_single'

    input:
    tuple val(clstr_num), val(clstr_order)
    path(clstr)

    output:
    path("${clstr_num}_sub.clstr"),emit:"subclstr"

    script:

    """

    clstr_merge.pl ${clstr_order} > ${clstr_num}_sub.clstr

    """

}
