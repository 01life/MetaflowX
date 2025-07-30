
process MULTICHECKM2 {
    
    tag "$id"
    
    label 'process_high'

    input:
    tuple val(id),path(bin_folders)
    path(checkm2_db)
    val(step_name)

    output:
    tuple val(id),path("${id}*_quality_report.tsv"), emit: "new_quality_report", optional: true
    path("${id}_checkm2_error.txt"), emit: "checkm2_error", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.checkm2_options ?: ""

    """


    # Create output directory
    mkdir -p ${id}

    # Loop through each input folder
    for bin_folder in ${bin_folders}; do
        folder_name=\$(basename \$bin_folder)
        output_dir=${id}/\${folder_name}

        checkm2 predict --threads ${task.cpus} --input \${bin_folder} --output-directory \${output_dir} --database_path ${checkm2_db} ${options} && cp \${output_dir}/quality_report.tsv ${id}-\${folder_name}_${step_name}_quality_report.tsv || {

cat <<-OUTLOG > \${output_dir}_checkm2_error.txt

==========Start at : `date` ==========
### Step ${task.process}
CheckM2 task for sample ${id}, folder \${folder_name} failed. Please check if bins are present in folder \${bin_folder}.
==========End at : `date` ==========

OUTLOG
        
        }
    done



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
        
    """
    
    stub:
    """
    mkdir -p ${id}
    touch ${id}1_quality_report.tsv
    """

}
