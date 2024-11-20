process REBINASSEMBLY {

	tag "$bin_id"

	label 'process_high'

    input:
    tuple val(bin_id),path(bin_fa),path(fq1),path(fq2)

    output:
	tuple val(bin_id),path("${bin_id}/${bin_id}_reassembly_contigs.fa"), emit:"rebin"

    when:
    task.ext.when == null || task.ext.when

    script:
    def options = params.hybridSPAdes_options ?: ""
	def maxmem = task.memory.toGiga()


    """
	mkdir bra_tmp
	spades.py -o bra_tmp \\
		-1 ${fq1} \\
		-2 ${fq2} \\
		--untrusted-contigs ${bin_fa} \\
		-t ${task.cpus} \\
		--memory $maxmem \\
		${options} || { rm -rf bra_tmp && exit 1; } 
	
	mkdir ${bin_id} 
	mv bra_tmp/contigs.fasta ${bin_id}/${bin_id}_reassembly_contigs.fa
	
	rm -rf bra_tmp

    """

}


