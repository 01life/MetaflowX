
process BIGMAP {

    tag "$id"

    label 'process_single','error_ignore'

    input:
    tuple val(id),path(reads),path(antismash_out)
    path(bigspace_db)

    output:
    path("${id}_BiG-MAP/csv-results/${id}_BiG-MAP_*.xls"),emit:"BGCfamily"
    path("${id}_BiG-MAP/csv-results/${id}_BiG-MAP_corecov.xls"),emit:"corecov"
    path("${id}_BiG-MAP/csv-results/*")

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_args = params.single_end ? "-U ${reads}" : "-I1 ${reads[0]} -I2 ${reads[1]}"
    def bigmap_family_options = params.bigmap_family_options ?: ""
    def bigmap_map_options = params.bigmap_map_options ?: "" 
    """

    mkdir all.bin_family
    BiG-MAP.family.py -D ${antismash_out} -b ${params.bigspace_path} -pf ${bigspace_db} -O all.bin_family/ -p ${task.cpus} ${bigmap_family_options}
    

    BiG-MAP.map.py ${reads_args} -O ${id}_BiG-MAP -F all.bin_family -th ${task.cpus} ${bigmap_map_options}

    simplify_antimash_ID_group.abundance.py ${id}_BiG-MAP/csv-results/BiG-MAP.map.results.ALL.csv ${id} ${id}_BiG-MAP/csv-results/${id}_BiG-MAP
    
    rm -rf ${id}_BiG-MAP/bowtie2* ${id}_BiG-MAP/bedtools-results ${id}_BiG-MAP/biom-results
    rm -rf all.bin_family
    
    """

    stub:
    """
    mkdir ${id}_BiG-MAP/csv-results
    touch ${id}_BiG-MAP/csv-results/${id}_BiG-MAP_corecov.xls
    
    """

}
