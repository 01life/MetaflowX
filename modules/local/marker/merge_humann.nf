
process MERGEHUMANN {

    label 'process_single'

    input:
    val(sample_number)
    val(finish_number)
    path("tables/*")
    path(samplesheet)


    output:
    path("${params.pipeline_prefix}*.xls"),emit: "humann_profile"
    path("humann_report"),emit: "humann_report"

    when:
    sample_number == finish_number

    script:
    def file_name = params.humann_version=="3" ? "pathcoverage" : "reactions"
    """

    # No sorting for single sample
    if [ ${sample_number} == 1 ]; then

        parallel "cp tables/*{1}.xls ${params.pipeline_prefix}_HUMAnN_{1}.xls" ::: genefamilies pathabundance ${file_name} 

        parallel "cp tables/*{1}_relab.xls ${params.pipeline_prefix}_HUMAnN_{1}_relab.xls ; humann_split_stratified_table -i ${params.pipeline_prefix}_HUMAnN_{1}_relab.xls -o . ; rm -f ${params.pipeline_prefix}_HUMAnN_{1}_relab.xls" :::  genefamilies eggnog go ko level4ec MetaCyc pfam

    else

        humann_join_tables -i tables -o ${params.pipeline_prefix}_HUMAnN_genefamilies_raw.txt --file_name genefamilies.xls
    col_reorder.pl ${params.pipeline_prefix}_HUMAnN_genefamilies_raw.txt 1 1 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_HUMAnN_genefamilies.xls
    
        parallel "humann_join_tables -i tables -o ${params.pipeline_prefix}_HUMAnN_{1}_raw.txt --file_name {1}; col_reorder.pl ${params.pipeline_prefix}_HUMAnN_{1}_raw.txt 1 1 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_HUMAnN_{1}.xls" ::: pathabundance ${file_name} 

        parallel "humann_join_tables -i tables -o ${params.pipeline_prefix}_HUMAnN_{1}_relab_raw.txt --file_name {1}_relab; col_reorder.pl ${params.pipeline_prefix}_HUMAnN_{1}_relab_raw.txt 1 1 <(sed 's/,/\\t/g' ${samplesheet}) 1 1 1 ${params.pipeline_prefix}_HUMAnN_{1}_relab.xls;humann_split_stratified_table -i ${params.pipeline_prefix}_HUMAnN_{1}_relab.xls -o . ; rm -f ${params.pipeline_prefix}_HUMAnN_{1}_relab.xls" :::  genefamilies eggnog go ko level4ec MetaCyc pfam
        
        # Delete result files before sorting.
        rm -rf ${params.pipeline_prefix}_HUMAnN_*_raw.txt

    fi
    
    sed -i '1s/[^\\t]*//' ${params.pipeline_prefix}_HUMAnN*relab*stratified.xls

    sed -i '1s/[^\\t]*//' ${params.pipeline_prefix}_HUMAnN_MetaCyc_relab*stratified.xls
    
    mkdir humann_report
    ln -s ../${params.pipeline_prefix}_HUMAnN_MetaCyc_relab_unstratified.xls humann_report/metacyc.xls
           

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS

    """

    stub:
    """
    touch ${params.pipeline_prefix}_HUMAnN_eggnog.xls
    touch ${params.pipeline_prefix}_HUMAnN_go.xls
    touch ${params.pipeline_prefix}_HUMAnN_ko.xls
    touch ${params.pipeline_prefix}_HUMAnN_level4ec.xls
    touch ${params.pipeline_prefix}_HUMAnN_pfam.xls
    touch ${params.pipeline_prefix}_HUMAnN_MetaCyc_relab_unstratified.xls        
    mkdir humann_report
    cp ${params.pipeline_prefix}_HUMAnN_MetaCyc_relab_unstratified.xls humann_report/metacyc.xls
    """ 

}
