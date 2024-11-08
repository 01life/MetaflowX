//
// SUBWORKFLOW: Read in samplesheet, validate and stage input files
//

include { READSVALIDATERENAME4PE as RAWREADSVALIDATERENAME4PE } from '../../modules/local/common/reads_validate_rename'
include { READSVALIDATERENAME4SE as RAWREADSVALIDATERENAME4SE } from '../../modules/local/common/reads_validate_rename'
include { READSVALIDATERENAME4PE as CLEANREADSVALIDATERENAME4PE } from '../../modules/local/common/reads_validate_rename'
include { READSVALIDATERENAME4SE as CLEANREADSVALIDATERENAME4SE } from '../../modules/local/common/reads_validate_rename'

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

def checkContent(LinkedHashMap row, String type) {
    
    //Determine if sequence data is empty.
    def id = row.id
    def contig = row.contig
    def fq1 = row[type.concat("_reads1")] ? file(row[type.concat("_reads1")]) : false
    def fq2 = row[type.concat("_reads2")] ? file(row[type.concat("_reads2")]) : false
    def se = row[type.concat("_se")] ? file(row[type.concat("_se")]) : false
    if(id == null) {exit 1, "Invalid input samplesheet: id can not be empty !"}
    
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

        if (!se) {  exit 1, "Invalid input samplesheet: ${type}_se can not be empty for ${id}!"  }
        
        return [ id, se, contig ]

    }else{

        if (!fq1) {  exit 1, "Invalid input samplesheet: ${type}_reads1 can not be empty for ${id}!"  }
        if (!fq2) {  exit 1, "Invalid input samplesheet: ${type}_reads2 can not be empty for ${id}!"  }
        
        return [ id, fq1, fq2, contig ]

    }

}

workflow INPUT_CHECK {
    take:
    type
    samplesheet     // file: /path/to/samplesheet.csv

    main:

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
        raw_reads_original = input_rows.map{ it ->  //[id, [raw_reads1, raw_reads2]] or [id, [raw_se]]
                        if (params.single_end)
                            return [ it[0], [it[1]] ]
                        else
                            return [ it[0], [it[1], it[2]] ]
                        }
                        
        // Raw reads validate and rename
        raw_reads = Channel.empty()
        if (params.single_end){
            RAWREADSVALIDATERENAME4SE(type, raw_reads_original)
            raw_reads = RAWREADSVALIDATERENAME4SE.out.reads
        }else{
            RAWREADSVALIDATERENAME4PE(type, raw_reads_original)
            raw_reads = RAWREADSVALIDATERENAME4PE.out.reads
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
        // Clean reads validate and rename
        clean_reads = Channel.empty()
        if(params.single_end){
            CLEANREADSVALIDATERENAME4SE(type, clean_reads_original)
            clean_reads = CLEANREADSVALIDATERENAME4SE.out.reads
        }else{
            CLEANREADSVALIDATERENAME4PE(type, clean_reads_original)
            clean_reads = CLEANREADSVALIDATERENAME4PE.out.reads
        }
        
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

