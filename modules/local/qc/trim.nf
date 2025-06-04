process TRIM {

    tag "$id"
    label 'process_single','error_ignore'

    input:
    tuple val(id), path(reads)
    path(host_db)
    path(adapters)
    path(phix_db)

    output:
    tuple val(id), path("${id}*clean*.fq.gz"), emit: "clean_reads"
    path("${id}.md5")
    tuple val(id), path("${id}_raw.txt"), path("${id}_trimmomatic.log"), path("${id}_rmhost.log"), emit: "stat_info"
    path("${id}_rmphix.log", optional: true)

    when:
    task.ext.when == null || task.ext.when

    script:
    def trim_options = params.trim_options ?: ""
    def ILLUMINACLIP = "ILLUMINACLIP:${adapters}:${params.trim_ILLUMINACLIP_options}"
    def bowtie2_options = params.qc_bowtie2_options ? "--bowtie2-options \"${params.qc_bowtie2_options}\"" : ""

    def trimCmd = ""
    def removeCmd = ""
    def cleanupCmd = ""

    def rmhostInput = params.single_end
        ? "-i1 ${id}_SE_trimmed.fq"
        : "-i1 ${id}_paired_1.fq -i2 ${id}_paired_2.fq"

    def phixInputFromHostClean = params.single_end
        ? "-i1 ${id}_host_rmhost_1.fq"
        : "-i1 ${id}_host_rmhost_1.fq -i2 ${id}_host_rmhost_2.fq"

    def fqToRemove = params.single_end
        ? "${id}_SE_trimmed.fq"
        : "${id}_paired_*.fq ${id}_unpaired_*.fq"

    def tempCleanFqToRemove = (params.host_db && params.phix_db) ? "${id}_host_rmhost_*.fq" : ""

    // Build trimming command
    if (params.single_end) {
        trimCmd = """
        fastq_stat.pl -s ${id} -a ${reads} -o ${id}_raw.txt

        trimmomatic SE -threads ${task.cpus} ${reads} ${id}_SE_trimmed.fq \\
            ${ILLUMINACLIP} ${trim_options} \\
            2> ${id}_trimmomatic.log
        """
    } else {
        trimCmd = """
        fastq_stat.pl -s ${id} -a ${id}_raw_1.fq.gz -b ${id}_raw_2.fq.gz -o ${id}_raw.txt

        trimmomatic PE -threads ${task.cpus} ${id}_raw_1.fq.gz ${id}_raw_2.fq.gz \\
            ${id}_paired_1.fq ${id}_unpaired_1.fq \\
            ${id}_paired_2.fq ${id}_unpaired_2.fq \\
            ${ILLUMINACLIP} ${trim_options} \\
            2> ${id}_trimmomatic.log
        """
    }

    // Host-only removal
    if (host_db && !phix_db) {
        removeCmd = """
        rmhost_with_bowtie2.py ${rmhostInput} -db ${host_db}/${params.host_db_index} \\
            -t ${task.cpus} -o ${id}_host ${bowtie2_options} 2>${id}_rmhost.log
        rename host_rmhost clean ${id}_host_rmhost*.fq
        """
    }

    // Phix-only removal
    if (phix_db && !host_db) {
        removeCmd = """
        rmhost_with_bowtie2.py ${rmhostInput} -db ${phix_db}/${params.phix_db_index} \\
            -t ${task.cpus} -o ${id}_phix ${bowtie2_options} 2>${id}_rmphix.log
        rename phix_rmhost clean ${id}_phix_rmhost*.fq
        """
    }

    // Both host and phix removal in sequence
    if (host_db && phix_db) {
        removeCmd = """
        rmhost_with_bowtie2.py ${rmhostInput} -db ${host_db}/${params.host_db_index} \\
            -t ${task.cpus} -o ${id}_host ${bowtie2_options} 2>${id}_rmhost.log

        rmhost_with_bowtie2.py ${phixInputFromHostClean} -db ${phix_db}/${params.phix_db_index} \\
            -t ${task.cpus} -o ${id}_phix ${bowtie2_options} 2>${id}_rmphix.log

        rename phix_rmhost clean ${id}_phix_rmhost*.fq
        """
    }

    // Final cleanup
    cleanupCmd = """
    pigz ${id}_clean*.fq
    rm -f ${fqToRemove} ${tempCleanFqToRemove}
    md5sum ${id}_clean*.fq.gz > ${id}.md5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """

    return """
    ${trimCmd}
    ${removeCmd}
    ${cleanupCmd}
    """
}
