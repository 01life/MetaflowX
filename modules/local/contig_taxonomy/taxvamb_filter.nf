
process TAXVAMBFILTER {
    
    tag "$id"
    
    label 'process_single'

    input:
    tuple val(id),path(bins)
    path(checkm2_db)

    output:
    tuple val(id),path("${id}/quality_report.tsv"), emit: "quality_report", optional: true
    tuple val(id),path("${id}_taxvamb_filter"), emit: "bins", optional: true
    path("${id}_checkm2_error.txt"), emit: "checkm2_error", optional: true
    tuple val(id),path("taxvamb.contigs2bin.tsv"), emit:"tsv", optional: true
    tuple val(id),path("${id}_TaxVamb_BinsContigs.tsv"), emit:"BinsContigs", optional: true



    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.checkm2_options ?: ""

    """

    checkm2 predict --threads ${task.cpus} --input ${bins} --output-directory ${id} --database_path ${checkm2_db} ${options} || {
    
    cat <<-OUTLOG > ${id}_checkm2_error.txt

    ==========Start at : `date` ==========
    ### Step ${task.process}
    CheckM2 task for sample ${id} failed. Please check if bins are present in folder ${bins}.
    ==========End at : `date` ==========

    OUTLOG
    
    }


    awk -F "\\t" '{if(NR==1)print \$0"\\t""QS"}{if(NR>1){sum=0;sum=\$2-5*\$3;print \$0"\\t"sum}}' ${id}/quality_report.tsv > ${id}_QS_quality_report.tsv
    
    awk -F "\\t" '{if ((\$2 > ${params.completeness})&&(\$3 < ${params.contamination})&&(\$NF > ${params.QS})) {print \$1".fa"}}' ${id}_QS_quality_report.tsv > filter.txt
    
    line_count=\$(wc -l < filter.txt)

    if [ \$line_count -gt 0 ]; then

        mkdir ${id}_taxvamb_filter
        for file in \$(cat filter.txt); do
            cp ${bins}/\$file ${id}_taxvamb_filter
        done

        Fasta_to_Contig2Bin.sh -i ${id}_taxvamb_filter -e fa > taxvamb.contigs2bin.tsv

        awk -F "\\t" '{print\$2"\\t"\$1"\\tTaxVamb"}' taxvamb.contigs2bin.tsv  > ${id}_TaxVamb_BinsContigs.tsv

    
    else
        
        cat <<-OUTLOG > ${id}_bin_filter.log
        
==========Start at : `date` ==========
### Step ${task.process}
Filtering bins based on bin completeness threshold of ${params.completeness} and contamination threshold of ${params.contamination}, yielded no results in sample ${id}.
==========End at : `date` ==========

OUTLOG

    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
        
    """

}
