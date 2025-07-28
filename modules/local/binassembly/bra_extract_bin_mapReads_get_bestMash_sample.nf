
process EXTRACTREADSMASHSAMPLE {

	tag "$bin_id"

	label 'process_single'

	input:
	tuple val(bin_id),val(top_sample_list),val(ref_fa),val(native_sample),path(reads1),path(reads2),path(bin_fa) //[binis,sample_txt(;),ref,nativesample,[reads1,2 list],bin_fa]


	output:
	tuple val(bin_id),path("BRA_${bin_id}/${bin_id}_bwa_mash_1.fq.gz"),path("BRA_${bin_id}/${bin_id}_bwa_mash_2.fq.gz"),emit:"map"
	tuple val(bin_id),path("BRA_${bin_id}/${bin_id}*Top_*_sample.txt"),emit:"binsample"

	when:
	task.ext.when == null || task.ext.when

	script:

	def options = params.extract_bin_reads_options ?: ""

	"""
	
	mkdir clean_fq
	mv ${reads1} clean_fq
	mv ${reads2} clean_fq
	

	bra_preprocess_bin_assembly.py -i ${bin_id} -a ${native_sample} -1 ${top_sample_list} -f clean_fq -b ${bin_fa} -r ${ref_fa} --threads ${task.cpus} ${options} BRA_${bin_id} 

	rm -rf BRA_${bin_id}/MappingFq


	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
		Python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python //; s/ .*\$//')
	END_VERSIONS

	"""
	stub:
	"""
	mkdir BRA_${bin_id}
	touch BRA_${bin_id}/${bin_id}_bwa_mash_1.fq.gz
	touch BRA_${bin_id}/${bin_id}_bwa_mash_2.fq.gz
	touch BRA_${bin_id}/${bin_id}_Top_10_sample.txt
	"""

}
