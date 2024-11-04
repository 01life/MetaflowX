//
// SUBWORKFLOW: Read in samplesheet, validate and stage input files
//

include { CLEANREADSRENAME } from '../../modules/local/common/clean_reads_rename'

//Check for duplicate sample IDs.
def checkID(file){

    def ID = file.readLines().drop(1).collect { it.split(",")[0] }
    
    // Total ID count.
    id_count = ID.size() 
    // Unique ID count.
    unique_id_count = ID.unique().size()

    if (id_count != unique_id_count) {
        throw new IllegalArgumentException("Duplicate sample ID found, please check input file: ${file}")           
    }

    return file

}

// Check if the file is a valid gz file.
def isGzFile(file) {  
    
    // Check if it is a file and its file extension.
    if(!(file.getName().endsWith(".gz")) || !file.isFile()) {
        return false
    }
    
    // The magic number for gzip files is 0x1F 0x8B.
    def magicBytes = new byte[2]  
    file.withInputStream { stream ->  
        stream.read(magicBytes)  
        return magicBytes[0] == (byte) 0x1F && magicBytes[1] == (byte) 0x8B  
    }  

}

// Check fastq files.
def checkFastq(filePath, String type, String suffix) {  
    // Fastq file path is not configured.
    if (!filePath) {  
        exit 1, "Invalid input samplesheet: ${type}_${suffix} can not be empty!"  
    // Check if the configured fastq files are valid.
    } else if (!isGzFile(filePath)) {  
        exit 1, "Error: ${filePath} is not a valid gz file, please check!"  
    }  
}

def checkContent(LinkedHashMap row, String type) {
    
    //Determine if sequence data is empty.
    def id = row.id
    def contig = row.contig
    def fq1 = row[type.concat("_reads1")] ? file(row[type.concat("_reads1")], checkIfExists: true) : false
    def fq2 = row[type.concat("_reads2")] ? file(row[type.concat("_reads2")], checkIfExists: true) : false
    def se = row[type.concat("_se")] ? file(row[type.concat("_se")], checkIfExists: true) : false
    if(id == null) {exit 1, "Invalid input samplesheet: id can not be empty !"}
    
    // SE data check for raw_se/clean_se.
    if(params.single_end){
        checkFastq(se, type, "se")
    // PE data check.
    }else{
        checkFastq(fq1, type, "reads1")
        checkFastq(fq2, type, "reads2")
    }
    
    //In modes 4/5, check the contig column.
    if( params.mode==4 || params.mode==5){
        if(contig != null){ 
            //Check if the file exists.
            if ( !file(contig).exists() ) {
                exit 1, "Contig file of sample ${id} does not exist !"
            }
        } else { exit 1, "Invalid input samplesheet: contig can not be empty !" }
    }
    
    //Check if sampleID is valid (letters, numbers, underscores).
    if ( !id.matches('^[a-zA-Z0-9_]*$') ) {
        exit 1, "Sample ID: ${id} contains invalid characters!"
    }
    

    if(params.single_end){
        return [ id, se, contig ]
    }else{
        //Check the size of the sequence files.
        def fq1size = fq1.size()
        def fq2size = fq2.size()
        if (Math.abs(fq1size - fq2size) > 2 * Math.min(fq1size, fq2size)) {
            exit 1, "The file sizes of reads1 and reads2 sequence files for sample ${id} differ significantly !"
        }
        return [ id, fq1, fq2, contig ]
    }

}

workflow INPUT_CHECK {
    take:
    type
    samplesheet     // file: /path/to/samplesheet.csv

    main:

    //If there were errors in the previous execution, clear the log to avoid re-running without issues and not generating a report.
    // def error_log = file("${params.outdir}/Stop_pipeline_error.log")
    // if(error_log.exists() && error_log.size()>0) {
    //     error_log.write("")
    // }

    //Check for duplicate IDs; exit immediately if duplicates are found.
    freads = checkID(samplesheet)

    //Parse and check inputs.
    input_rows = Channel.from ( freads )
        .splitCsv ( header:true, sep:',' )
        .map { row ->
            if(row.size()>=2){
                checkContent( row, type )
            }else{
                exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects at least 2 !"
            }
        }

    //Pass sequence data to the corresponding channel based on reads type.
    if(type == "raw"){
        raw_reads = input_rows.map{ it ->  //[id, [raw_reads1, raw_reads2]] or [id, [raw_se]]
                        if (params.single_end)
                            return [ it[0], [it[1]] ]
                        else
                            return [ it[0], [it[1], it[2]] ]
                        }
        clean_reads = Channel.empty()
        contig = Channel.empty()
    }else{
        raw_reads = Channel.empty()
        clean_reads_original = input_rows.map{ it -> //[id, [clean_reads1, clean_reads2]] or [id, clean_se]
            if (params.single_end)
                return [ it[0], [it[1]] ]
            else
                return [ it[0], [it[1], it[2]] ]
            }                    

        // For input clean reads, check if the filename meets the requirements.
        CLEANREADSRENAME(clean_reads_original)
        clean_reads = CLEANREADSRENAME.out.reads


        // Parse contig data in modes 4/5.
        if (params.mode in [ 4, 5 ]){
            contig = input_rows.map{ it -> [it[0], it[-1]] } // [id, contig]
        }else{
            contig = Channel.empty()
        }

    }


    emit:
    raw_reads           // channel: [ val(id), [ raw reads1, raw reads2 ] ]
    clean_reads           // channel: [ val(id), [ clean reads1, clean reads2 ] ]
    contig           // channel: [ val(id), [ contig ] ]

}

