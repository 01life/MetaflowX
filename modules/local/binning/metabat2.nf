
process METABAT2 {
    
    tag "$id"
    
    label 'process_single','error_ignore'
    
    input:
    tuple val(id),path(contigs),path(bin_depth)

    output:
    tuple val(id),path("bins",type:'dir'),emit:"bins"
    tuple val(id),path("metabat.contigs2bin.tsv"),emit:"tsv"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.metabat2_options ?: ""
    """

    metabat2 -t ${task.cpus} -i ${contigs} -a ${bin_depth} -o bins/${id} ${options}

    finish=0
    if [ -d "bins" ] ; then      
        finish=\$((ls -1 "./bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        
        cd bins/
        for file in *.fa; do mv "\$file" "metabat2_\${file}"; done
        cd ..

        Fasta_to_Contig2Bin.sh -i bins -e fa > metabat.contigs2bin.tsv
    
    else
        touch metabat.contigs2bin.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS

    """
}
