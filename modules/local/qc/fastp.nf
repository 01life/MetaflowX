process FASTP {

    tag "$id"
    label 'process_single', 'error_ignore'

    input:
    tuple val(id), path(reads) 
    path(host_db)
    path(adapters)
    path(phix_db)

    output:
    tuple val(id), path("${id}*clean*.fq.gz"), emit: "clean_reads"
    path("${id}_fastp.json"), emit: "fastp_json"
    path("${id}_fastp.html")
    path("${id}.md5")
    tuple val(id), path("${id}_fastp.json"), path("${id}_rmhost.log"), emit: "stat_info"
    path("${id}_rmphix.log", optional: true)

    when:
    task.ext.when == null || task.ext.when

    script:
    def fastpOptions = params.fastp_options ?: ""
    def bowtie2Options = params.qc_bowtie2_options ? "--bowtie2-options \"${params.qc_bowtie2_options}\"" : ""

    def fastpCmd = ""
    def removeCmd = ""

    def rmhostInput = params.single_end
        ? "-i1 ${id}.fq"
        : "-i1 ${id}_paired_1.fq -i2 ${id}_paired_2.fq"

    def phixInputFromHostClean = params.single_end
        ? "-i1 ${id}_host_rmhost_1.fq"
        : "-i1 ${id}_host_rmhost_1.fq -i2 ${id}_host_rmhost_2.fq"

    def fqToRemove = params.single_end
        ? "${id}.fq"
        : "${id}_paired_*.fq"

    def tempCleanFqToRemove = (params.host_db && params.phix_db) ? "" : "${id}_host_rmhost_*.fq"

    // Build fastp command
    fastpCmd = params.single_end
        ? """
        fastp \\
            --in1 ${reads} \\
            --out1 ${id}.fq \\
            --thread ${task.cpus} \\
            --adapter_fasta ${adapters} \\
            ${fastpOptions}
        """
        : """
        fastp \\
            --in1 ${id}_raw_1.fq.gz \\
            --in2 ${id}_raw_2.fq.gz \\
            --out1 ${id}_paired_1.fq \\
            --out2 ${id}_paired_2.fq \\
            --thread ${task.cpus} \\
            --adapter_fasta ${adapters} \\
            ${fastpOptions}
        """

    // Host-only removal
    if (host_db && !phix_db) {
        removeCmd = """
        rmhost_with_bowtie2.py ${rmhostInput} -db ${host_db}/${params.host_db_index} \\
            -t ${task.cpus} -o ${id}_host ${bowtie2Options} 2>${id}_rmhost.log
        rename host_rmhost clean ${id}_host_rmhost*.fq
        """
    }

    // Phix-only removal
    if (phix_db && !host_db) {
        removeCmd = """
        rmhost_with_bowtie2.py ${rmhostInput} -db ${phix_db}/${params.phix_db_index} \\
            -t ${task.cpus} -o ${id}_phix ${bowtie2Options} 2>${id}_rmphix.log
        rename phix_rmhost clean ${id}_phix_rmhost*.fq
        """
    }

    // Both host and phix removal in sequence
    if (host_db && phix_db) {
        removeCmd = """
        rmhost_with_bowtie2.py ${rmhostInput} -db ${host_db}/${params.host_db_index} \\
            -t ${task.cpus} -o ${id}_host ${bowtie2Options} 2>${id}_rmhost.log

        rmhost_with_bowtie2.py ${phixInputFromHostClean} -db ${phix_db}/${params.phix_db_index} \\
            -t ${task.cpus} -o ${id}_phix ${bowtie2Options} 2>${id}_rmphix.log
        rename phix_rmhost clean ${id}_phix_rmhost*.fq
        """
    }

    // Cleanup and metadata generation
    cleanupCmd = """
    pigz ${id}_clean*.fq
    rm -f ${fqToRemove} ${tempCleanFqToRemove}
    md5sum ${id}_clean*.fq.gz > ${id}.md5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """


    return """
    ${fastpCmd}
    rename fastp ${id}_fastp fastp.*

    ${removeCmd}

    ${cleanupCmd}
    """
}
