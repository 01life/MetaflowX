 
include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'
include { GTDB } from '../../modules/local/binning/gtdb'
include { MERGEGTDB } from '../../modules/local/binning/merge_gtdb'
include { GTDBWARNING } from '../../modules/local/binning/gtdb_warning'

workflow BINTAXONOMY {
    take:
    bins_folder         // RENAMEBIN.out.folder
    bins_rename_map     // RENAMEBIN.out.name_mapping
    qs_quality_report   //POSTBINNING.out.qs_quality_report
    sample_bin_map

    main:
    
    /*
    * Verify the essential parameters for running this module
    */
    bin_taxon_essential_db = [params.gtdbtk_db, params.mash_db, params.gtdb_archaeal_metadata, params.gtdb_bacterial_metadata]

    if(!checkEssentialParams(bin_taxon_essential_db)) { exit 1, "The required parameters to execute the Bin Taxonomy module are:\n --gtdbtk_db\n --mash_db\n --gtdb_archaeal_metadata\n --gtdb_bacterial_metadata" }

    ch_gtdbtk_db = params2Channel(params.gtdbtk_db)
    ch_mash_db = params2Channel(params.mash_db)
    ch_gtdb_archaeal_metadata = params2Channel(params.gtdb_archaeal_metadata)
    ch_gtdb_bacterial_metadata = params2Channel(params.gtdb_bacterial_metadata)


    ch_gtdb_input = bins_folder.flatten()
                            .map{ chunk ->
                                [ chunk.baseName,chunk ]                                
                                }
    GTDB (ch_gtdb_input, ch_gtdbtk_db, ch_mash_db, ch_gtdb_archaeal_metadata, ch_gtdb_bacterial_metadata)

    ch_gtdb2ncbi = GTDB.out.gtdb2ncbi.collectFile(name:"gtdb_taxonomy_to_ncbi.xls", keepHeader:true)
    
    ch_gtdbtk_summary = GTDB.out.gtdbtk_summary.collectFile(name:"merged_gtdbtk_summary.xls", keepHeader:true)
    
    summary_linecount = ch_gtdbtk_summary.countLines()

    //Execute if the number of rows in the merged result file is greater than 1.
    MERGEGTDB(summary_linecount, ch_gtdbtk_summary, ch_gtdb2ncbi, qs_quality_report, bins_rename_map, sample_bin_map)
    gtdb_report = MERGEGTDB.out.gtdb_report
    gtdb_summary = MERGEGTDB.out.gtdb_summary
    bin_QS_taxonomy = MERGEGTDB.out.bin_QS_taxonomy

    //Execute if the number of rows in the merged result file is less than 2.
    GTDBWARNING(summary_linecount)

    gtdb_report = gtdb_report.mix(GTDBWARNING.out.log).collect()

    emit:
    gtdb_report
    gtdb_summary
    bin_QS_taxonomy

}