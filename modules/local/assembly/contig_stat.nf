
process CONTIGSTAT {

    label 'process_low'
    
    input:
    val(sample_number)
    val(finish_number)
    path(all_contig)

    output:
    path("contig.path.list")
    path("contig_report"),emit:"contig_report"
    path("all_contig_info.txt"),emit:"contig_info"

    when:
    sample_number == finish_number

    script:
    """

    cat ${all_contig} > all_contig
    contig_stat.py all_contig all_contig_info.txt

    rm -rf all_contig

    mkdir contig_report
    cp all_contig_info.txt contig_report/contigstat.txt
    
    readlink -f *fa |awk -F "/" '{print\$NF"\\t"\$0}' |sed 's/_contigs.fa\\t/,/g' > contig.path.list

    """

}
