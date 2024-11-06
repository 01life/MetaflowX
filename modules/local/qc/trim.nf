
process TRIM {
    
    tag "$id"

    label 'process_single','error_ignore'

    input:
    tuple val(id),path(reads)
    path(host_db)
    path(adapters)

    output:
    path("${id}.md5")
    tuple val(id),path("${id}*.fq.gz"),emit:"clean_reads"
    tuple val(id),path("${id}_raw.txt"),path("${id}_rmhost.log"),emit:"stat_info"

    when:
    task.ext.when == null || task.ext.when

    script:
    def trim_options = params.trim_options ?: ""
    def ILLUMINACLIP = "ILLUMINACLIP:${adapters}:${params.trim_ILLUMINACLIP_options}"
    def bowtie2_options = params.qc_bowtie2_options ? "--bowtie2-options \"${params.qc_bowtie2_options}\"" : ""
    
    if(params.single_end){
        
        """        
        #Raw reads information statistics.
        fastq_stat.pl -s ${id} -a ${reads} -o ${id}_raw.txt

        trimmomatic SE -threads ${task.cpus} ${reads} ${id}_SE_paired.fq ${ILLUMINACLIP} ${trim_options}

        rmhost_with_bowtie2.py -i1 ${id}_SE_paired.fq -db ${host_db}/${params.host_db_index} -t ${task.cpus} -o ${id} ${bowtie2_options} 2>${id}_rmhost.log

        rename rmhost clean ${id}_rmhost*.fq
        pigz ${id}_clean*.fq
        rm -rf ${id}*.fq

        md5sum ${id}_clean*.fq.gz > ${id}.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            java: \$(java --version | sed 's/java //g')
        END_VERSIONS
            
        """
    }else{

        """
        #Raw reads information statistics.
        fastq_stat.pl -s ${id} -a ${id}_1.fq.gz -b ${id}_2.fq.gz -o ${id}_raw.txt

        trimmomatic PE -threads ${task.cpus} ${id}_raw_1.fq.gz ${id}_raw_2.fq.gz ${id}_paired_1.fq ${id}_unpaired_1.fq ${id}_paired_2.fq ${id}_unpaired_2.fq ${ILLUMINACLIP} ${trim_options}

        rmhost_with_bowtie2.py -i1 ${id}_paired_1.fq -i2 ${id}_paired_2.fq -db ${host_db}/${params.host_db_index} -t ${task.cpus} -o ${id} ${bowtie2_options} 2>${id}_rmhost.log

        rename rmhost clean ${id}_rmhost*.fq
        pigz ${id}_clean*.fq
        rm -rf ${id}*.fq

        md5sum ${id}_clean_*.fq.gz > ${id}.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimmomatic: \$(trimmomatic -version)
        END_VERSIONS

        """
    }

}
