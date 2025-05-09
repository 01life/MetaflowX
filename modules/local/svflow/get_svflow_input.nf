
process GETSVFLOWINPUT {

    label 'process_low'

    input:
    path(all_profile)
    path(profile_csv)

    output:
    path("${params.pipeline_prefix}_svflow_input.csv"), optional: true
    path("${params.pipeline_prefix}_svflow_input_warning.log"), emit: "log", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    # Iterate through input files, filter data.
    while IFS=, read -r col1 col2 col3; do
    # Check if the second column file exists.
    if [[ -f "\$col2" ]]; then
      # Write data to the output file.
      echo "\$col1,\$col2,\$col3" >> ${params.pipeline_prefix}_svflow_input.csv
    fi
    done < ${profile_csv}

    if [ ! -f ${params.pipeline_prefix}_svflow_input.csv ] ; then
      echo "The input csv for SVflow is empty, please check your data !" > ${params.pipeline_prefix}_svflow_input_warning.log
    fi

    """
}
