
process MULTICOV {

    tag "$binID"

    label 'process_single'

    input:
    tuple val(binID),path(covList)


    output:
    tuple val(binID), path("${binID}_mulit_cov.txt"),emit:"mulitCov"
    

    when:
    task.ext.when == null || task.ext.when

    script:
   
    def cov_list = covList instanceof List ? covList.join(" ") : "$covList"

    """
    bra_merge_converage.py  ${cov_list} > ${binID}_mulit_cov.txt
    """

    stub:
    """
    touch ${binID}_mulit_cov.txt
    """
}
