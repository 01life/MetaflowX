
process GETSAMPLEBINMAP {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(bins,stageAs:'tmp/*')

    output:
    path("${id}_bin_mapping.txt"), emit:"mapping"

    when:
    task.ext.when == null || task.ext.when

    script:

    // def bin_list = bins instanceof List ? bins.join("/*.fa ") : "$bins/*.fa"

    """

    ls tmp/*/*.fa |awk -F "/" '{print"${id}\\t"\$NF}' > ${id}_bin_mapping.txt

    """

    stub:
    """
    echo -e "${id}\tbin1" > ${id}_bin_mapping.txt
    """
}
