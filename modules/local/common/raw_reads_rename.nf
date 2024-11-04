
process RAWREADSRENAME {

    tag "$id"

    label 'process_low'

    input:
    tuple val(id),path(reads)

    output:
    tuple val(id),path("${id}/*.fq.gz"),emit:"reads"

    when:
    task.ext.when == null || task.ext.when

    script:
    if(params.single_end){
        
        """
        mkdir ${id}
        ln -sf \$(readlink -f ${reads}) ${id}/${id}.fq.gz
    
        """

    }else{

        """
        mkdir ${id}
        ln -sf \$(readlink -f ${reads[0]}) ${id}/${id}_1.fq.gz
        ln -sf \$(readlink -f ${reads[1]}) ${id}/${id}_2.fq.gz

        """
    }

}
