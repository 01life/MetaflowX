
process CATCONTIG  {
    
    tag "$id"

    label 'process_high'
    
    input:
    tuple val(id),path(contigs),path(proteins_fasta)
    path(cat_gtdb_db)
    path(cat_pack)

    output:
    tuple val(id), path("CAT_${id}_run.summary.txt"),emit:"cat_report"
    tuple val(id), path("${id}_CAT_run.contig2classification_filter_withoutScores.txt"),emit:"cat_taxonomy"
    tuple val(id), path("${id}_CAT_run.contig2classification_filter_withoutScores.txt"),emit:"cat_taxonomy_score"

    when:
    task.ext.when == null || task.ext.when

    script:
    def cat_contig_options = params.cat_contig_options ?: ""

    def  prodigal_faa = params.mode in [5, 0]  ? " -p ${proteins_fasta} " : ""


    """
    export PATH=${cat_pack}:\$PATH

    CAT_pack contigs -c ${contigs} -d ${cat_gtdb_db}/db -t ${cat_gtdb_db}/tax -n ${task.cpus} --out_prefix ${id}_CAT_run ${prodigal_faa} ${cat_contig_options}

    filter_cat_annotation.py -i ${id}_CAT_run.contig2classification.txt -o ${id}_CAT_run.contig2classification_filter_withScores.txt -O ${id}_CAT_run.contig2classification_filter_withoutScores.txt



 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$(echo \$(CAT_pack --version 2>&1) | awk '{print\$3}')
    END_VERSIONS

    """

    stub:
    """
    touch CAT_${id}_run.summary.txt

    """

}
