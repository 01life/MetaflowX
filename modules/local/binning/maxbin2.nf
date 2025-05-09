process MAXBIN2 {
   
    tag "$id"
   
    label 'process_single','error_ignore'
    
    input:
    tuple val(id),path(contigs),path(reads)

    output:
    tuple val(id),path("bins",type:'dir'), emit:"bins", optional: true
    tuple val(id),path("maxbin.contigs2bin.tsv"), emit:"tsv", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_args = params.single_end ? "-reads ${reads}" : "-reads ${reads[0]} -reads2 ${reads[1]}"
    def options = params.maxbin2_options ?: ""
    
    """
    mkdir bins
    run_MaxBin.pl -contig ${contigs} -out bins/${id} ${reads_args} -thread ${task.cpus} ${options} || echo "MAXBIN2 task for sample ${id} failed ......" > ${id}.log

    finish=0
    if [ -d "bins" ] && [ ! -e "${id}.log" ] ; then      
        finish=\$((ls -1 "./bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        
        cd bins
        for file in *.fasta; do mv "\$file" "maxbin2_\${file%.fasta}.fa"; done
        cd ..

        Fasta_to_Contig2Bin.sh -i bins -e fa > maxbin.contigs2bin.tsv
    
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maxbin2: \$( run_MaxBin.pl -v | head -n 1 | sed 's/MaxBin //' )
    END_VERSIONS

    """
}
