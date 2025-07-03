process METADECODER {

    tag "$id"

    label 'process_single','error_ignore'

    input:
    tuple val(id), path(contigs), path(sorted_bam)

    output:
    tuple val(id), path("${id}_metadecoder_bins/", type: 'dir'), emit: "bins", optional: true
    tuple val(id), path("MetaDecoder.contigs2bin.tsv"), emit: "tsv", optional: true
    tuple val(id), path("${id}_MetaDecoder_BinsContigs.tsv"), emit: "BinsContigs", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.metadecoder_options ?: ""

    """
    metadecoder coverage --threads ${task.cpus} -b ${sorted_bam} -o ${id}.metadecoder.coverage

    metadecoder seed --threads ${task.cpus} -f ${contigs} -o ${id}.metadecoder.seed

    metadecoder cluster ${options} -f ${contigs} -c ${id}.metadecoder.coverage -s ${id}.metadecoder.seed -o ${id}.metadecoder || echo "MetaDecoder task for sample ${id} failed ......" > ${id}.log

    finish=0
    if [ ! -e "${id}.log" ]; then
        files=( ${id}.metadecoder.*.fasta )
        if [ -e "\${files[0]}" ]; then
            finish=\${#files[@]}
        fi
    fi

    if [ \$finish -gt 0 ]; then

        mkdir ${id}_metadecoder_bins
        mv ${id}.metadecoder.*.fasta ./${id}_metadecoder_bins

        cd ${id}_metadecoder_bins

        i=1
        for file in *; do
            new_name="Metadecoder_${id}_bin.\$i.fa"
            mv "\$file" "\$new_name"
            i=\$((i+1))
        done

        cd ../

        Fasta_to_Contig2Bin.sh -i ${id}_metadecoder_bins -e fa > MetaDecoder.contigs2bin.tsv

        awk -F "\\t" '{print \$2"\\t"\$1"\\tMetaDecoder"}' MetaDecoder.contigs2bin.tsv > ${id}_MetaDecoder_BinsContigs.tsv

    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaDecoder: \$( metadecoder --version )
    END_VERSIONS
    """
}
