
process MERGEGENEINFO {

    label 'process_low'

    input:
    val(sample_number)
    val(finish_number)
    path(all_stat_info)

    output:
    path("*.multi.cdhit.task"),emit:"multi_task" , optional: true
    path("single.cdhit.task"),emit:"single_task" , optional: true
    path("geneset_sampleStat_report"),emit:"report"


    when:
    sample_number == finish_number

    script:
    
    """

    mkdir geneset_sampleStat_report
    cp ${all_stat_info} geneset_sampleStat_report/genesetSampleStat.txt

    check_cdhit_task_num.py -i ${all_stat_info} -s ${params.cdhit_split_run_threshold} -c ${params.cdhit_geneset_chunk_size} 


    """

}
