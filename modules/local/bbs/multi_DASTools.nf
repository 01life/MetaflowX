
process MULTIDASTOOL {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(contigs),path(faa),val(label),path(contig2bin)

    output:
    path("${id}_${label}_allBins.eval.txt"), emit: "eval", optional: true
    path("${id}_DAS_Tool_error.txt"), emit: "das_bins_error", optional: true
    path("${id}_${label}_contigs2bin.tsv"), emit: "tsv", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def bin_list = contig2bin instanceof List ? contig2bin.join(",") : "$contig2bin"
    def dastool_options = params.dastool_options ?: "" 
    """
    
    DAS_Tool -i ${bin_list} -c ${contigs} -o ${label} -p ${faa} -t ${task.cpus} ${dastool_options} || {
            
    cat <<-OUTLOG > ${id}_${label}_DAS_Tool_error.txt

    ==========Start at : `date` ==========
    ### Step ${task.process}
    DASTool task for sample ${id} failed. It is spossible that no bins provided, please check contig2bin files: ${bin_list}. Alternatively, no bins found with scores greater than 0.5(default). This can occur if the input data is of poor quality or if the parameters used to run DASTool are not optimal.
    ==========End at : `date` ==========

    OUTLOG

    }
    
    if ls "${label}_DASTool_bins/"*.fa 1> /dev/null 2>&1; then
        Fasta_to_Contig2Bin.sh -i ${label}_DASTool_bins/ -e fa > ${id}_${label}_contigs2bin.tsv
        sed 's/^/${id}\\t${label}\\t/g' ${label}_allBins.eval  > ${id}_${label}_allBins.eval.txt
    fi
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dastool: \$( DAS_Tool --version 2>&1 | grep "DAS Tool" | sed 's/DAS Tool //' )
    END_VERSIONS
    """

}
