
process READSVALIDATERENAME4SE {

    tag "$id"

    label 'process_low'

    input:
    val(tag)
    tuple val(id),path(reads)

    output:
    tuple val(id),path("${id}/*.fq.gz"),emit:"reads"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    header=\$(hexdump -n 2 -e '1/1 "%02x"' "${reads}")  
    if [ "\$header" == "1f8b" ]; then  
        echo "${reads} is a valid gz file."  
    else  
        echo "${reads} is not a valid gz file or it is corrupted."  
        exit 1  
    fi  

    mkdir ${id}
    ln -sf \$(readlink -f ${reads}) ${id}/${id}_${tag}.fq.gz

    """

}


process READSVALIDATERENAME4PE {

    tag "$id"

    label 'process_low'

    input:
    val(tag)
    tuple val(id),path(reads)

    output:
    tuple val(id),path("${id}/*.fq.gz"),emit:"reads"

    when:
    task.ext.when == null || task.ext.when

    script:
    def fq1size = reads[0].size()
    def fq2size = reads[1].size()
    def diff = Math.abs(fq1size - fq2size)
    def min = 2 * Math.min(fq1size, fq2size)

    """
    # Check the size of the sequence files
    if [ ${diff} -gt ${min} ]; then
        echo "The file sizes of reads1 and reads2 sequence files for sample ${id} differ significantly !"
        exit 1 
    fi

    # Check if the file is a valid gz file
    header1=\$(hexdump -n 2 -e '1/1 "%02x"' "${reads[0]}")  
    if [ "\$header1" == "1f8b" ]; then  
        echo "${reads[0]} is a valid gz file."  
    else  
        echo "${reads[0]} is not a valid gz file or it is corrupted."  
        exit 1  
    fi  

    header2=\$(hexdump -n 2 -e '1/1 "%02x"' "${reads[1]}")  
    if [ "\$header2" == "1f8b" ]; then  
        echo "${reads[1]} is a valid gz file."  
    else  
        echo "${reads[1]} is not a valid gz file or it is corrupted."  
        exit 1  
    fi  

    mkdir ${id}
    ln -sf \$(readlink -f ${reads[0]}) ${id}/${id}_${tag}_1.fq.gz
    ln -sf \$(readlink -f ${reads[1]}) ${id}/${id}_${tag}_2.fq.gz

    """

}
