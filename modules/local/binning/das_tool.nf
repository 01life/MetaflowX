
process DASTOOL {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(contigs),path(faa),path(bins)

    output:
    tuple val(id),path("${id}_DASTool_bins",type:'dir'), emit: "das_bins", optional: true
    path("${id}_DAS_Tool_error.txt"), emit: "das_bins_error", optional: true
    path("${id}*.tsv"), optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def bin_list = bins instanceof List ? bins.join(",") : "$bins"
    def dastool_options = params.dastool_options ?: "" 
    """
    
    DAS_Tool -i ${bin_list} -c ${contigs} -o ${id} -p ${faa} -t ${task.cpus} ${dastool_options} || {
    
    cat <<-OUTLOG > ${id}_DAS_Tool_error.txt

    ==========Start at : `date` ==========
    ### Step ${task.process}
    DASTool task for sample ${id} failed. It is spossible that no bins provided, please check contig2bin files: ${bin_list}. Alternatively, no bins found with scores greater than 0.5(default). This can occur if the input data is of poor quality or if the parameters used to run DASTool are not optimal.
    ==========End at : `date` ==========

    OUTLOG
    
    }

    #Process output during normal execution.
    if ls "${id}_DASTool_bins/"*.fa 1> /dev/null 2>&1; then
        mv ${id}_allBins.eval ${id}_allBins_eval.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dastool: \$( DAS_Tool --version 2>&1 | grep "DAS Tool" | sed 's/DAS Tool //' )
    END_VERSIONS
    """

}
