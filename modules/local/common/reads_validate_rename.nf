process READSVALIDATERENAME {
    
    tag "$id"

    label 'process_low'

    input:
    val(tag)
    tuple val(id),path(reads)

    output:
    tuple val(id),path("${id}/*.fq.gz"),emit:"reads",optional:true
    path("${id}_error.log"),emit:"error",optional:true
    path("${id}_check.log"),emit:"info",optional:true

    when:
    task.ext.when == null || task.ext.when

    script:
    def fq_count = reads instanceof List ? reads.size() : 1
    """
    function check_gz_file {
        local file="\$1"
        
        if [ ! -e "\$(readlink -f "\$file")" ]; then
            echo "ERROR: \$file: No such file." >> ${id}_error.log
            return 1
        else
            echo "PASS: \$file exist." >> ${id}_check.log
        fi

        local header=\$(hexdump -n 2 -e '1/1 "%02x"' "\$file")
        if [ "\$header" != "1f8b" ]; then
            echo "Sample ${id}: \$file is not a valid gz file or it is corrupted." >> ${id}_error.log
            return 1
        else 
            echo "Sample ${id}: \$file is a valid gz file." >> ${id}_check.log
            return 0
        fi
    }
    
    function validate_file_sizes() {  
        local file1=\$1  
        local file2=\$2  
        fq1size=\$(stat -c%s \$(readlink -f "\$file1"))  
        fq2size=\$(stat -c%s \$(readlink -f "\$file2"))  
        diff=\$((fq1size > fq2size ? fq1size - fq2size : fq2size - fq1size))  
        min=\$((2 * (fq1size < fq2size ? fq1size : fq2size)))  

        if [ \$diff -gt \$min ]; then  
            echo "There's a notable size difference between the reads1 and reads2 sequence files for sample ${id}, please check!" >> ${id}_error.log
            return 1
        else
            echo "There's a tiny difference in file size between your reads1 and reads2 sequence files for sample ${id}. Nothing to worry about, but just so you know!" >> ${id}_check.log
            return 0
        fi  
    } 

    mkdir ${id}

    # Check and rename for single-end or first read of paired-end
    if [ ${fq_count} -eq 1 ]; then

        (check_gz_file "${reads}" && ln -sf \$(readlink -f "${reads}") ${id}/${id}_${tag}.fq.gz) || echo "${reads} error"

    else
        (check_gz_file "${reads[0]}" &&
        check_gz_file "${reads[1]}" && 
        validate_file_sizes "${reads[0]}" "${reads[1]}" &&
        ln -sf \$(readlink -f "${reads[0]}") ${id}/${id}_${tag}_1.fq.gz &&
        ln -sf \$(readlink -f "${reads[1]}") ${id}/${id}_${tag}_2.fq.gz) ||
        echo "${reads} error"
        
    fi

    """
}