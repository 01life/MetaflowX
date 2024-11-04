
process METABINNER {
    
    tag "$id"
    
    label 'process_medium','error_ignore'

    input:
    tuple val(id),path(contigs),path(bin_depth)

    output:
    tuple val(id),path("${id}/metabinner_bins/",type:'dir'),emit:"bins"
    tuple val(id),path("metabinner.contigs2bin.tsv"),emit:"tsv"


    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.metabinner_options ?: ""
    def length =  params.min_contig_len - 1

    """
    ln -s ${contigs} contigs.fa
    
    python ${params.metabinner_path}/scripts/gen_kmer.py contigs.fa ${length} 4

    cut -f 1,4 ${bin_depth} > coverage_profile.tsv

    run_metabinner.sh \
        -p ${params.metabinner_path} \
        -a \$PWD/${contigs} \
        -t ${task.cpus} \
        -d \$PWD/coverage_profile.tsv \
        -k \$PWD/contigs_kmer_4_f${length}.csv \
        -o \$PWD/${id} ${options}

    recover_binning_pro.py -t \$PWD/${id}/metabinner_res/metabinner_result.tsv -f ${contigs} -o \$PWD/${id}/metabinner_bins/ 
    
    finish=0
    if [ -d "${id}/metabinner_bins" ] ; then      
        finish=\$((ls -1 "./${id}/metabinner_bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then

        cd ${id}/metabinner_bins/
        for file in *.fa; do mv "\$file" "metabinner_${id}_\${file}"; done
        cd ../..

        Fasta_to_Contig2Bin.sh -i ${id}/metabinner_bins/ -e fa > metabinner.contigs2bin.tsv

        #sed -i 's%\\tbin%\\t${id}%g' metabinner.contigs2bin.tsv
    
    else
        touch metabinner.contigs2bin.tsv
    fi


    rm -rf ${id}/metabinner_res/intermediate_result ${id}/metabinner_res/ensemble_res ${id}/metabinner_res/unitem_profile

    """

}
