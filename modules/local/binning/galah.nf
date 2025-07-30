
process GALAH {

    label 'process_high'

    input:
    path(final_bins,stageAs:'bins/*')
    path(checkm2_report,stageAs:'report/*' )

    output:
    path("Galah_cluster/*.fa"),emit:"genomes"

    when:
    task.ext.when == null || task.ext.when

    script:

    def galah_options = params.galah_options ?: "" 
    
    """
    awk 'FNR==1 && NR!=1 {next;}{print}'  report/*tsv  > all.checkm2-quality-report.tsv

    galah cluster --genome-fasta-directory bins -x fa --output-representative-fasta-directory Galah_cluster --checkm2-quality-report all.checkm2-quality-report.tsv ${galah_options}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Galah: \$(echo \$(galah --version 2>&1) | sed 's/galah //g' )
    END_VERSIONS 

    """

    stub:
    """
    mkdir -p Galah_cluster
    touch Galah_cluster/representative.fa
    touch all.checkm2-quality-report.tsv
    """

}
