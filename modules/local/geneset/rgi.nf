
process RGI {

    label 'process_single'

    input:
    path(pep)
    path(abundance)
    path(CARD_db)

    output:
    path("${params.pipeline_prefix}_geneset_function_CARD*.xls")
    path("${params.pipeline_prefix}_geneset_function_CARD_annotation.xls"),emit:"rgi_anno"
    path("${params.pipeline_prefix}_geneset_function_RGI_annotation.xls"),emit:"rgi"
    path("rgi_report"),emit:"rgi_info"

    when:
    task.ext.when == null || task.ext.when

    script:
    def rgi_load_options = params.rgi_load_options ?: ""
    def rgi_main_options = params.rgi_main_options ?: ""
    
    """
        
    rgi clean --local

    #Load the database.
    rgi load \
    --card_json ${CARD_db}/card.json \
    --card_annotation ${CARD_db}/card_database_v3.2.7.fasta \
    --card_annotation_all_models ${CARD_db}/card_database_v3.2.7_all.fasta \
    --wildcard_annotation ${CARD_db}/wildcard_database_v3.2.7.fasta \
    --wildcard_annotation_all_models ${CARD_db}/wildcard_database_v3.2.7_all.fasta \
    --wildcard_index ${CARD_db}/wildcard/index-for-model-sequences.txt \
    --amr_kmers ${CARD_db}/wildcard/all_amr_61mers.txt \
    --kmer_database ${CARD_db}/wildcard/61_kmer_db.json  \
    --local \
    ${rgi_load_options}
        
    #Protein annotation.
    sed 's/*//g' ${pep} > new.protein.fa
    
    rgi main --input_sequence new.protein.fa \
    --output_file ./RGI --num_threads ${task.cpus} \
    ${rgi_main_options}
            
    cut -f 1,9 RGI.txt |sed 's/ /_/g' |sed 's/,/_/g'|tail -n +2 > ${params.pipeline_prefix}_geneset_function_CARD_annotation.xls

    cp RGI.txt ${params.pipeline_prefix}_geneset_function_RGI_annotation.xls

    
    functional_profile.pl -k ${params.pipeline_prefix}_geneset_function_CARD_annotation.xls -f ${abundance} -o ./ -pr ${params.pipeline_prefix}_geneset_function_CARD_abundance.xls

    #Clean up the loaded database after execution.
    rgi clean --local

    mkdir rgi_report
    cp ${params.pipeline_prefix}_geneset_function_CARD_abundance.xls rgi_report/CARD.xls

    #echo CARD \$PWD/${params.pipeline_prefix}_geneset_function_CARD_abundance.xls > rgi.report.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$( rgi main --version )
    END_VERSIONS

    """

}
