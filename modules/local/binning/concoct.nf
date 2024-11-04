
process CONCOCT {
    
    tag "$id"
    
    label 'process_single','error_ignore'

    input:
    tuple val(id),path(contigs),path(sorted_bam),path(sorted_bam_csi)

    output:
    // tuple val(id),path("${id}_concoct_output/fasta_bins/*.fa"),emit:"fa"
    // path("${id}_concoct_output/",type:'dir')
    tuple val(id),path("${id}_concoct_output/fasta_bins",type:'dir'),emit:"bins"
    tuple val(id),path("concoct.contigs2bin.tsv"),emit:"tsv"


    when:
    task.ext.when == null || task.ext.when

    script:
    def concoct_options = params.concoct_options ?: ""
    
    """
    cut_up_fasta.py ${contigs} -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
    concoct_coverage_table.py contigs_10K.bed ${sorted_bam} > coverage_table.tsv
    concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b ${id}_concoct_output/ -t ${task.cpus} ${concoct_options}

    merge_cutup_clustering.py ${id}_concoct_output/clustering_gt1000.csv > ${id}_concoct_output/clustering_merged.csv
    mkdir ${id}_concoct_output/fasta_bins
    extract_fasta_bins.py ${contigs} ${id}_concoct_output/clustering_merged.csv --output_path ${id}_concoct_output/fasta_bins

    finish=0
    if [ -d "${id}_concoct_output/fasta_bins" ] ; then      
        finish=\$((ls -1 "./${id}_concoct_output/fasta_bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then

        cd ${id}_concoct_output/fasta_bins
        i=1; for file in *; do new_name="concoct_${id}_bin.\$i.fa"; mv "\$file" "\$new_name"; i=\$((i+1)); done
        cd ../..

        Fasta_to_Contig2Bin.sh -i ${id}_concoct_output/fasta_bins -e fa > concoct.contigs2bin.tsv

        #perl -pe \"s/,/\t${id}./g;\" ${id}_concoct_output/clustering_merged.csv > concoct.contigs2bin.tsv
        #sed -i '1d' concoct.contigs2bin.tsv

    else
        touch concoct.contigs2bin.tsv
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2>&1) | sed 's/concoct //g' )
    END_VERSIONS

    """

}
