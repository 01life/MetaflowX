
process HUMANNV4 {
    
    tag "$id"
    
    label 'process_high'

    input:
    tuple val(id),path(reads),path(profile)
    path(humann_chocophlan_db)
    path(humann_protein_db)
    path(humann_map_db)

    output:
    path("${id}/${id}*"),emit:"profile"
    tuple val(id),path("${id}/${id}_2_genefamilies.xls"),emit:"gf_profile"
    tuple val(id),path("${id}/${id}_ko_relab.xls"),emit:"ko_profile",optional:true

    when:
    task.ext.when == null || task.ext.when

    script:
    def pre_command = params.single_end ? "zcat ${reads} > ${id}.fq" : "zcat ${reads[0]} ${reads[1]} > ${id}.fq"
    def options1 = params.humann_options ?: ""
    def options2 = params.humann_regroup_table_options ?: ""
    """

    ${pre_command}

    humann -i ${id}.fq -o ${id} --threads ${task.cpus} --o-log ${id}/${id}.log --nucleotide-database ${humann_chocophlan_db} --taxonomic-profile ${profile} --protein-database ${humann_protein_db} ${options1}

    #humann_renorm_table -i ${id}/${id}_2_genefamilies.tsv -u relab -o ${id}/${id}_genefamilies_relab.tsv

    #parallel "humann_regroup_table ${options2} -i ${id}/${id}_genefamilies_relab.tsv -c ${humann_map_db}/map_{1}_uniref90.txt.gz -o ${id}/${id}_{1}_relab.tsv" ::: go pfam level4ec ko eggnog

    #humann_renorm_table -i ${id}/${id}_4_pathabundance.tsv -u relab -o ${id}/${id}_MetaCyc_relab.tsv
    
    rename .tsv .xls ${id}/*

    rm -rf ${id}.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$( humann --version | sed 's/humann //' )
    END_VERSIONS
        
    """

    stub:
    """

    mkdir ${id}
    touch ${id}/${id}_genefamilies_relab.xls
    touch ${id}/${id}_ko_relab.xls

    """
}
