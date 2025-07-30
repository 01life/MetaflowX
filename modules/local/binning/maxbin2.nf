process MAXBIN2 {
   
    tag "$id"
   
    label 'process_single','error_ignore'
    
    input:
    tuple val(id),path(contigs),path(reads)

    output:
    tuple val(id),path("maxbin2bins",type:'dir'), emit:"bins", optional: true
    tuple val(id),path("maxbin.contigs2bin.tsv"), emit:"tsv", optional: true
    tuple val(id),path("${id}_MaxBin2_BinsContigs.tsv"), emit:"BinsContigs", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_args = params.single_end ? "-reads ${reads}" : "-reads ${reads[0]} -reads2 ${reads[1]}"
    def options = params.maxbin2_options ?: ""
    
    """
    mkdir maxbin2bins
    run_MaxBin.pl -contig ${contigs} -out maxbin2bins/${id} ${reads_args} -thread ${task.cpus} ${options} || echo "MAXBIN2 task for sample ${id} failed ......" > ${id}.log

    finish=0
    if [ -d "maxbin2bins" ] && [ ! -e "${id}.log" ] ; then      
        finish=\$((ls -1 "./maxbin2bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        
        cd maxbin2bins
        #for file in *.fasta; do mv "\$file" "maxbin2_\${file%.fasta}.fa"; done
        i=1; for file in *.fasta; do new_name="maxbin2_${id}_bin.\$i.fa"; mv "\$file" "\$new_name"; i=\$((i+1)); done

        cd ..

        Fasta_to_Contig2Bin.sh -i maxbin2bins -e fa > maxbin.contigs2bin.tsv

        awk -F "\\t" '{print\$2"\\t"\$1"\\tMaxBin2"}' maxbin.contigs2bin.tsv  > ${id}_MaxBin2_BinsContigs.tsv
        
    
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maxbin2: \$( run_MaxBin.pl -v | head -n 1 | sed 's/MaxBin //' )
    END_VERSIONS

    """
    stub:
    """
    mkdir maxbin2bins
    touch maxbin2bins/${id}.log
    touch maxbin.contigs2bin.tsv
    touch ${id}_MaxBin2_BinsContigs.tsv
    """
}
