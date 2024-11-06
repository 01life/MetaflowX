
process FASTP {
    
    tag "$id"

    label 'process_single','error_ignore'

    input:
    tuple val(id),path(reads)
    path(host_db)
    path(adapters)

    output:
    tuple val(id),path("${id}*.fq.gz"),emit:"clean_reads"
    path("${id}_fastp.json"),emit:"fastp_json"
    path("${id}_fastp.html")
    path("${id}.md5")
    tuple val(id),path("${id}_fastp.json"),path("${id}_rmhost.log"),emit:"stat_info"

    when:
    task.ext.when == null || task.ext.when

    script:
    def fastp_options = params.fastp_options ?: ""
    def bowtie2_options = params.qc_bowtie2_options ? "--bowtie2-options \"${params.qc_bowtie2_options}\"" : ""

    if(params.single_end){
        """

        fastp \\
        --in1 ${reads} \\
        --out1 ${id}.fq \\
        --thread ${task.cpus} \\
        --adapter_fasta ${adapters} \\
        ${fastp_options}

        rename fastp ${id}_fastp fastp.* 

        rmhost_with_bowtie2.py -i1 ${id}.fq -db ${host_db}/${params.host_db_index} -t ${task.cpus} -o ${id} ${bowtie2_options} 2>${id}_rmhost.log

        rename rmhost clean ${id}_rmhost*.fq
        pigz ${id}_clean*.fq
        rm -rf ${id}.fq

        md5sum ${id}_clean*.fq.gz > ${id}.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
    
        """

    }else{

        """

        fastp \\
        --in1 ${id}_raw_1.fq.gz \\
        --out1 ${id}_paired_1.fq \\
        --in2 ${id}_raw_2.fq.gz \\
        --out2 ${id}_paired_2.fq \\
        --thread ${task.cpus} \\
        --adapter_fasta ${adapters} \\
        ${fastp_options}

        rename fastp ${id}_fastp fastp.* 

        rmhost_with_bowtie2.py -i1 ${id}_paired_1.fq -i2 ${id}_paired_2.fq -db ${host_db}/${params.host_db_index} -t ${task.cpus} -o ${id} ${bowtie2_options} 2>${id}_rmhost.log

        rename rmhost clean ${id}_rmhost*.fq
        pigz ${id}_clean*.fq
        rm -rf ${id}_paired*.fq

        md5sum ${id}_clean_*.fq.gz > ${id}.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        
        """
    }

}
