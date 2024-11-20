
process GETBINMAPREADS {

	tag "$bin_id"

	label 'process_single'

	publishDir(
        path: "${params.outdir}/07.BinReassembly",
        pattern: "BRA_${bin_id}/${bin_id}*",
        mode: "copy",
        failOnError: true
    )

    input:
    path(clean_fq)
    tuple val(bin_id),path(bin_fa),path(ref_fa),val(ha_fq_id),val(fq_id_list)


    output:
    tuple val(bin_id),path("BRA_${bin_id}/${bin_id}_bwa_mash_1.fq.gz"),path("BRA_${bin_id}/${bin_id}_bwa_mash_2.fq.gz"),emit:"map"

    when:
    task.ext.when == null || task.ext.when

    script:

    def options = params.preprocess_bin_assembly_options ?: ""

    """
    
    mkdir clean_fq
    mv ${clean_fq} clean_fq

    preprocess_bin_assembly.py -i ${bin_id} -a ${ha_fq_id} -1 ${fq_id_list} -f clean_fq -b ${bin_fa} -r ${ref_fa} --threads ${task.cpus} ${options} BRA_${bin_id} 

	rm -rf BRA_${bin_id}/MappingFq


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python //; s/ .*\$//')
    END_VERSIONS

    """

}
