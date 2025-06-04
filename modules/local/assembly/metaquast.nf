process METAQUAST {
    tag "$id"

    label 'process_medium'

    input:
    tuple val(id), path(contigs)

    output:
    path "${id}.tar.gz"                     , emit: qc
    path "${id}_transposed_report.tsv"      , emit: report
    // path "versions.yml"                  , emit: versions

    script:

    def metaquastOptions = params.metaquast_options ?: ""


    """
    metaquast.py --threads "${task.cpus}" -l "${id}" "${contigs}" -o "${id}" ${metaquastOptions}
    cp ${id}/transposed_report.tsv ${id}_transposed_report.tsv

    tar -zcf ${id}.tar.gz ${id}
    rm -rf ${id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        metaquast: \$(metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//")
    END_VERSIONS
    """
}
