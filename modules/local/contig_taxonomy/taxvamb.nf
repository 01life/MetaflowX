
process TAXVAMB { 
    
    tag "$id"

    label 'process_medium'
    
    input:
    tuple val(id),path(contigs),path(depth),path(taxonomy)

    output:
    // tuple val(id), path("${id}/"),emit:"taxvamb" , optional: true
    tuple val(id), path("${id}/results_taxometer.tsv"),emit:"taxometer" , optional: true
    // tuple val(id), path("${id}/vaevae_clusters_unsplit.tsv"),emit:"vaevae" , optional: true


    tuple val(id),path("${id}/TaxVamb",type:'dir'), emit:"bins", optional: true
    tuple val(id),path("taxvamb.contigs2bin.tsv"), emit:"tsv", optional: true
    tuple val(id),path("${id}_TaxVamb_BinsContigs.tsv"), emit:"BinsContigs", optional: true


    when:
    task.ext.when == null || task.ext.when

    script:
    def taxvamb_options = params.taxvamb_options ?: ""
    """
    cut -f 1,2  ${taxonomy} >  ${id}.taxonomy.tmp

    vamb bin taxvamb --outdir ${id} --fasta ${contigs} --abundance_tsv ${depth} --taxonomy ${id}.taxonomy.tmp -p ${task.cpus}   ${taxvamb_options} || echo "vamb bin taxvamb task for sample ${id} failed ......" > ${id}.log


    finish=0
    if [ -d "${id}" ] && [ ! -e "${id}.log" ]; then      
        finish=\$((ls -1 "${id}/bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        cd ${id}/bins
            i=1; for file in *; do new_name="taxvamb_${id}_bin.\$i.fa"; mv "\$file" "\$new_name"; i=\$((i+1)); done

        cd ..
            Fasta_to_Contig2Bin.sh -i bins -e fa > taxvamb.contigs2bin.tsv

            awk -F "\\t" '{print\$2"\\t"\$1"\\tTaxVamb"}' taxvamb.contigs2bin.tsv  > ${id}_TaxVamb_BinsContigs.tsv

        mv bins TaxVamb
        mv taxvamb.contigs2bin.tsv ../
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Vamb: \$( vamb --version )
    END_VERSIONS

    """

    stub:
    """
    mkdir ${id}-taxometer

    """

}
