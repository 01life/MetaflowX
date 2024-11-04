
process MERGEHUMANNEXPAND {

    tag "$name"

    label 'process_single'

    input:
    tuple val(name),path(abundance)
    path(samplesheet)

    output:
    path("${params.pipeline_prefix}*.xls")
    path("${name}.csv"),emit: "humann_profile"

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    mkdir tables
    mv ${abundance} tables

    humann_join_tables -i tables -o ${params.pipeline_prefix}_HUMAnN_${name}_raw.txt --file_name ${name}_relab

    col_reorder.pl ${params.pipeline_prefix}_HUMAnN_${name}_raw.txt 1 1 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_HUMAnN_${name}_relab.xls
    
    humann_split_stratified_table -i  ${params.pipeline_prefix}_HUMAnN_${name}_relab.xls -o ./
    
    rm -f ${params.pipeline_prefix}_HUMAnN_{1}_relab.xls ${params.pipeline_prefix}_HUMAnN_${name}_raw.txt

    echo ${name},${ch_output}/102.HUMAnN/${params.pipeline_prefix}_HUMAnN_${name}_relab_unstratified.xls,DAs > ${name}.csv

    """

}
