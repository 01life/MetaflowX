
process MULTIDASTOOL {
    
    tag "$id"
    
    label 'process_high'

    input:
    tuple val(id), path(dastoolinput),path(protein),path(contig)

    output:
    tuple val(id),path("${id}*_allBins.eval.txt"), emit: "eval", optional: true
    path("${id}*_DAS_Tool_error.txt"), emit: "das_bins_error", optional: true
    tuple val(id),path("${id}*_contigs2bin.tsv"), emit: "tsv", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def dastool_options = params.dastool_options ?: "" 
    """
    run_dastool_parallel.py --input ${dastoolinput} --cpus ${task.cpus}  --dastool_opts "${dastool_options}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dastool: \$( DAS_Tool --version 2>&1 | grep "DAS Tool" | sed 's/DAS Tool //' )
    END_VERSIONS
    """
    stub:
    """
    touch "${id}_allBins.eval.txt"
    touch "${id}_contigs2bin.tsv"
    """ 

}