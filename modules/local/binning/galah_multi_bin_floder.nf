
process GALAHMULTIBIN {

    tag "$id"

    label 'process_medium'

    input:
    tuple val(id),path(final_bins,stageAs:'multi_bins/*'),path(checkm2_report,stageAs:'report/*' )

    output:
    tuple val(id),path("${id}_Galah_cluster/"),emit:"genomes"
    tuple val(id),path("${id}_all.checkm2-quality-report_galah_selected.tsv"),emit:"report"

    when:
    task.ext.when == null || task.ext.when

    script:

    def galah_options = params.galah_options ?: "" 
    
    """
    mkdir -p galah_inputs
    
    cp multi_bins/*/*fa  ./galah_inputs/

    awk 'FNR==1 && NR!=1 {next;}{print}'  report/*tsv  > ${id}_all.checkm2-quality-report.tsv

    fa_count=\$(ls galah_inputs/*.fa | wc -l)

    if [ "\$fa_count" -le 1 ]; then
        echo "[INFO] Only one fasta file detected, skip galah cluster."
        mkdir -p ${id}_Galah_cluster
        single_fa=\$(ls galah_inputs/*.fa | head -n 1)
        cp "\$single_fa" ${id}_Galah_cluster/${id}_representative.fa

        cp ${id}_all.checkm2-quality-report.tsv ${id}_all.checkm2-quality-report_galah_selected.tsv

    else
        galah cluster --genome-fasta-directory galah_inputs -x fa --output-representative-fasta-directory ${id}_Galah_cluster --checkm2-quality-report ${id}_all.checkm2-quality-report.tsv ${galah_options}

        get_galah_result_qs.py -i ${id}_Galah_cluster -c ${id}_all.checkm2-quality-report.tsv -o ${id}_all.checkm2-quality-report_galah_selected.tsv
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Galah: \$(echo \$(galah --version 2>&1) | sed 's/galah //g' )
    END_VERSIONS 

    """
}
