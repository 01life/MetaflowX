
process TAXOMETER2BIN { 
    
    tag "$id"

    label 'process_low'
    
    input:
    tuple val(id),path(contigs),path(taxonomy)

    output:
    tuple val(id),path("tomomaBins",type:'dir'), emit:"bins", optional: true
    tuple val(id),path("tomoma.contigs2bin.tsv"), emit:"tsv", optional: true



    when:
    task.ext.when == null || task.ext.when

    script:

    """
    taxometer2Bin.py  -i ${taxonomy} -f ${contigs} -o tomomaBins -p tomoma_${id} || echo "TAXOMETER2BIN task for sample ${id} failed ......" > ${id}.log

    finish=0
    if [ -d "tomomaBins" ] && [ ! -e "${id}.log" ]; then      
        finish=\$((ls -1 "./tomomaBins") | wc -l)
    fi

    if [ \$finish -gt 0 ]; then
        
        Fasta_to_Contig2Bin.sh -i tomomaBins -e fa > tomoma.contigs2bin.tsv 
    
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$( python --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir tomomaBins
    touch  tomoma.contigs2bin.tsv 
    """

}
