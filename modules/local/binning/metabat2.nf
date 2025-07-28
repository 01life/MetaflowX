
process METABAT2 {
    
    tag "$id"
    
    label 'process_single','error_ignore'
    
    input:
    tuple val(id),path(contigs),path(bin_depth)

    output:
    tuple val(id),path("metabat2bins",type:'dir'), emit:"bins", optional: true
    tuple val(id),path("metabat.contigs2bin.tsv"), emit:"tsv", optional: true
    tuple val(id),path("${id}_MetaBAT2_BinsContigs.tsv"), emit:"BinsContigs", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.metabat2_options ?: ""
    """

    metabat2 -t ${task.cpus} -i ${contigs} -a ${bin_depth} -o metabat2bins/${id} ${options} || echo "METABAT2 task for sample ${id} failed ......" > ${id}.log

    finish=0
    if [ -d "metabat2bins" ] && [ ! -e "${id}.log" ]; then      
        finish=\$((ls -1 "./metabat2bins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        
        cd metabat2bins/
        #for file in *.fa; do mv "\$file" "metabat2_\${file}"; done
        i=1; for file in *; do new_name="metabat2_${id}_bin.\$i.fa"; mv "\$file" "\$new_name"; i=\$((i+1)); done

        cd ..

        Fasta_to_Contig2Bin.sh -i metabat2bins -e fa > metabat.contigs2bin.tsv

        awk -F "\\t" '{print\$2"\\t"\$1"\\tMetaBAT2"}' metabat.contigs2bin.tsv  > ${id}_MetaBAT2_BinsContigs.tsv

    
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS

    """
    stub:
    """
    mkdir -p metabat2bins
    touch metabat.contigs2bin.tsv
    touch ${id}_MetaBAT2_BinsContigs.tsv
    """
}
