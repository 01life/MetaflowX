
process DREP {

    label 'process_high'

    input:
    path(final_bins)

    output:
    path("drep_output/dereplicated_genomes/*.fa"),emit:"genomes"
    path("${params.pipeline_prefix}_dRep_cluster.xls")

    when:
    task.ext.when == null || task.ext.when

    script:
    def drep_split_thread = params.drep_split_thread ?: task.cpus
    def drep_split_mem = params.drep_split_mem ?: task.memory.toGiga()
    def options1 = params.drep_options ?: ""
    def options2 = params.drep_options ? "--drep_options=\"${params.drep_options}\"" : ""
    """
    facount=\$(ls -1 \$PWD/*.fa | wc -l)

    if [ \$facount -gt ${params.drep_split_run_threshold} ] && [ ${task.executor} != "local" ]; then

        drep_path=\$(which dRep)
        dRep_para.py drep_output -d \$drep_path -g ${final_bins} -s ${params.drep_bin_chunk_size} -S ${params.drep_split_run_threshold} -t ${drep_split_thread} -T ${task.executor} -m ${drep_split_mem} -p ${params.Account} -q ${task.queue} ${options2} 
        mv drep_output/dRep_cluster.txt ${params.pipeline_prefix}_dRep_cluster.xls

    else
        dRep dereplicate drep_output -g ${final_bins} -p ${task.cpus} ${options1}
        
        get_dRep_cluster.py ./drep_output/data_tables/Wdb.csv ./drep_output/data_tables/Cdb.csv ${params.pipeline_prefix}_dRep_cluster.xls

        #mv  drep_output/data_tables/Cdb.csv ${params.pipeline_prefix}_dRep_cluster.xls

    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dRep: \$( dRep 2>&1 | sed -n '2p' )
    END_VERSIONS 

    """

    stub:
    """
    mkdir -p drep_output/dereplicated_genomes
    touch drep_output/dereplicated_genomes/genome1.fa
    touch ${params.pipeline_prefix}_dRep_cluster.xls
    """

}
