process COMBINEBINNER {

    tag "$id"

    label 'process_single'

    input:
    tuple val(id), path(contigs2tsv), path(depth), path(contig)

    output:
    // Output the refined bin table (tsv), FASTA folder, and contig-to-bin mapping file
    tuple val(id), path("*__*tsv"), emit: "combine_bin_tsv", optional: true
    tuple val(id), path("*:*", type: 'dir'), emit: "combine_bin_fa", optional: true
    tuple val(id), path("*allcontigs2bin.txt"), emit: "allcontigs2bin", optional: true

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    set -e

    echo "Running comebine_bins_refine.py for sample ${id}"
    
    comebine_bins_refine.py -info ${depth} -bins ${contigs2tsv} -t ${task.cpus} -i ${id} -c ${contig}
    status=\$?

    # If the script failed or no valid output was produced, write an error log
    if [[ \$status -ne 0 || ! -s "${id}_allcontigs2bin.txt" ]]; then
    cat <<OUTLOG > ${id}_PBO_error.txt
    ========== Start at : \$(date) ==========
    ### Step: NFCORE_METASSEMBLY:METASSEMBLY:BINNING:BINNER:COMBINEBINNER
    The PBO process for combining multiple binning results for sample '${id}' has failed.

    Possible reasons:
    - No valid binning input files were provided. Please check: metabat.contigs2bin.tsv concoct.contigs2bin.tsv semibin2.contigs2bin.tsv metabinner.contigs2bin.tsv
    - No bin sets passed the internal quality check or threshold criteria.

    ========== End at : \$(date) ==========
    OUTLOG
    fi

    """

    stub:
    """
    mkdir -p MetaDecoder:semibin2
    mkdir -p semibin2:metabinner
    touch ${id}_allcontigs2bin.txt
    touch ${id}_MetaDecoder__semibin2.contigs2bin.tsv
    touch ${id}_semibin2__metabinner.contigs2bin.tsv
    """
    
}
