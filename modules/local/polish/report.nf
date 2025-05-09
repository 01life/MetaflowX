
process REPORT {

    label 'process_medium'

    input:
    path(result)
    path(log)
    path(topic)
    path(order)
    path(template)
    path(images)

    output:
    path("*.html"),emit:"html"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    
    if [ -n "\$(find . -maxdepth 1  -name 'coverm*' -print -quit)" ]; then
        echo "The directory exists."

        mkdir COVERM_REPORT
        cp -f coverm*/*xls COVERM_REPORT
        #rm -rf ./coverm*

        ls ./COVERM_REPORT/* |sed  -e ':x;N;s/\\n/,/;b x'  |sed 's/^/binabundance /g' > all.report

    fi

    ls  *_report/* |grep -v coverm | awk -F "/" '{print\$NF}' |awk -F "." '{print\$1}' > p0
    ls ./*_report/* |grep -v coverm > p1 
    paste -d " " p0 p1 >> all.report
    rm -rf p0 p1

    report_main_V20240509.py -r all.report -t ${topic} -i ${order} -c ${template} -p ${params.pipeline_prefix} -m ${images}

    if [ "${params.remove_temp_sam_bam}" == "true" ]; then
        find ${workDir} -name "*.sam" -type f -delete
        find ${workDir} -name "*.bam" -type f -delete
        find ${workDir} -name "*.bt2" -type f -delete
    fi

    """
}
