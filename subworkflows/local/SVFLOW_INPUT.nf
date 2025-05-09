
include { GETGCANDDEPTH } from '../../modules/local/svflow/get_contigs_gc_depth'
include { BUILDTREE } from '../../modules/local/svflow/build_tree'
include { COVERMRELABUNST } from '../../modules/local/svflow/coverm_rel_abun_st'
include { GETBINSTAXON } from '../../modules/local/svflow/get_bins_taxon'
include { GETBINSINFO } from '../../modules/local/svflow/get_final_bins_info'
include { MERGETREE } from '../../modules/local/svflow/merge_tree'
include { GETSVFLOWINPUT } from '../../modules/local/svflow/get_svflow_input'
include { PIPELINEWARNING } from '../../modules/local/common/pipeline_warning'


workflow SVFLOWINPUT {

    take:
    outdir
    contig_info     //ch_contig_info
    marker_profile  //RAPID_TAXONOMIC_PROFILING.out.ch_profile
    gene_func_abun  //GENEPROFILE.out.func_abun
    bins_list        //RENAMEBIN.out.list
    bins_folder     //RENAMEBIN.out.folder
    bins_count      //RENAMEBIN.out.count
    rawbin_quality  //POSTBINNING.out.qs_quality_report
    bins_info       //RENAMEBIN.out.bins_info
    rename_map      //RENAMEBIN.out.name_map
    gtdb_summary    //MERGEGTDB.out.gtdb_summary
    depth_list      //BOWTIE2BIN.out.list
    coverm_rel_abun //BINABUNDANCE.out.totalRelativeAbun


    main:

    ch_all_profile = Channel.empty()

    /*
    * marker + geneset result
    */
    ch_all_profile = ch_all_profile.mix(marker_profile, gene_func_abun)


    /*
    * Bins result analysis profile.
    */
    
    // gc depth
    // ch_all_depth = depth_list.collectFile(name:"all.depth.list")
    GETGCANDDEPTH(contig_info, depth_list.collect(), bins_list)
  
    ch_gtdb_input = bins_folder.flatten()
        .map{ chunk ->
            [ chunk.baseName,chunk ]                                
            }

    // GTDB tree construction; execute tree building task only if the total number of bins is less than or equal to gtdb_bin_chunk_size (default 500).
    BUILDTREE(bins_count, ch_gtdb_input) 

  
    // Normalize coverm results.
    COVERMRELABUNST(coverm_rel_abun)

    // Obtain bins classification information.
    GETBINSTAXON(gtdb_summary, COVERMRELABUNST.out.st_rel_abun)
    
    // Obtain bins information.
    GETBINSINFO(GETBINSTAXON.out.taxon, rawbin_quality, bins_info, rename_map)

    // Merge tree information.
    MERGETREE(BUILDTREE.out.tree, GETBINSTAXON.out.taxon, GETBINSINFO.out.bins_info)

    ch_all_profile = ch_all_profile.mix(GETGCANDDEPTH.out.info, COVERMRELABUNST.out.st_rel_abun, GETBINSTAXON.out.profile, GETBINSINFO.out.bins_info, MERGETREE.out.tree)



    /*
    * Merge all profile information.
    */
    ch_path_map = Channel.fromPath("${params.path_map}").splitCsv()

    ch_profile = Channel.fromPath("${params.profile_template}").splitCsv().groupTuple()

    ch_data = ch_path_map.join(ch_profile).transpose()
            .map{ it -> 
                return it[2] +","+ outdir + it[1] + "${params.pipeline_prefix}" + it[3] +","+ it[-1]
            }
            .collectFile(name:"profile.csv", newLine:true)

    GETSVFLOWINPUT(ch_all_profile.collect(), ch_data)

    PIPELINEWARNING("getSVflowInput", GETSVFLOWINPUT.out.log)

}