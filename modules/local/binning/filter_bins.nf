
process FILTERBINS {
    
    tag "$id"
    
    label 'process_low'

    input:
    tuple val(id),path(bins),path(quality_report)

    output:
    path("${id}/*.fa"),emit:"filter_bins", optional: true
    path("${id}_QS_quality_report.tsv"), emit: "QS_quality_report"
    path("${id}_bin_filter.log"), emit:"filtered_log", optional: true
    tuple val(id), path("${id}",type:'dir'), emit: "bins", optional: true
    tuple val(id), path("${id}_filtered_quality_report.tsv"), emit: "qs", optional: true
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """

    awk -F "\\t" '{if(NR==1)print \$0"\\t""QS"}{if(NR>1){sum=0;sum=\$2-5*\$3;print \$0"\\t"sum}}' ${quality_report} > ${id}_QS_quality_report.tsv
    

    awk -F "\\t" 'NR==1 {print \$0; next} {qs=\$2 - 5*\$3; if ((\$2 > '${params.completeness}') && (\$3 < '${params.contamination}') && (qs > '${params.QS}')) print \$0}' ${quality_report} > ${id}_filtered_quality_report.tsv


    awk -F "\\t" '{if ((\$2 > ${params.completeness})&&(\$3 < ${params.contamination})&&(\$NF > ${params.QS})) {print \$1".fa"}}' ${id}_QS_quality_report.tsv > filter.txt

    
    line_count=\$(wc -l < filter.txt)

    if [ \$line_count -gt 0 ]; then

        mkdir ${id}
        for file in \$(cat filter.txt); do
            cp ${bins}/\$file ${id}
        done
    
    else
        
    cat <<-OUTLOG > ${id}_bin_filter.log
            
    ==========Start at : `date` ==========
    ### Step ${task.process}
    Filtering bins based on bin completeness threshold of ${params.completeness} and contamination threshold of ${params.contamination}, yielded no results in sample ${id}.
    ==========End at : `date` ==========

    OUTLOG

    fi

    
    """
}
