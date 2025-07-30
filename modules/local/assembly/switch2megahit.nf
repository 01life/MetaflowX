
process SWITCH2MEGAHIT {

    tag "$id"

    label 'process_low'

    input:
    tuple val(id),path(contigs)

    output:
    path("${id}_assembly_warning.log"),emit:"log"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    cat <<-ASSEMBLYLOG > ${id}_assembly_warning.log
    
    ==========Start at : `date` ==========
    ### Step ${task.process}
    Sample ID: ${id} failed to assemble using metaSPAdes, so the assembly was completed using the MEGAHIT tool instead.
    ==========End at : `date` ==========

    ASSEMBLYLOG

    """
    
    stub:
    """
    touch ${id}_assembly_warning.log
    """

}
