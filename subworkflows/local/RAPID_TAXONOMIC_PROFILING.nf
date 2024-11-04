//
// metaphlan humann kraken
//

include { params2Channel } from '../../modules/local/common/utils'
include { checkEssentialParams } from '../../modules/local/common/utils'
include { METAPHLANV40 } from '../../modules/local/marker/metaphlan_v4_0'
include { METAPHLANV41 } from '../../modules/local/marker/metaphlan_v4_1'
include { MERGEMPA } from '../../modules/local/marker/merge_mpa'
include { SBG2GTDB } from '../../modules/local/marker/sbg2gtdb'
include { MPAEXTRAABUN } from '../../modules/local/marker/mpa_extra_abun'
include { MPAEXTRAABUN as MPAEXTRAABUN4HUMANN } from '../../modules/local/marker/mpa_extra_abun'
include { MERGEMPAEXTRAABUN } from '../../modules/local/marker/merge_mpa_extra_abun'
include { KRAKEN2 } from '../../modules/local/marker/kraken2'
include { MERGEKRAKEN2 } from '../../modules/local/marker/merge_kraken2'
include { HUMANNV3 } from '../../modules/local/marker/humann_v3'
include { HUMANNV4 } from '../../modules/local/marker/humann_v4'
include { MERGEHUMANN } from '../../modules/local/marker/merge_humann'
include { HUMANNEXPAND } from '../../modules/local/marker/humann_expand_DB'
include { MERGEHUMANNEXPAND } from '../../modules/local/marker/merge_humann_expand'
include { PIPELINEERROR as PIPELINEERROR_MPA} from '../../modules/local/common/pipeline_error'
include { PIPELINEERROR as PIPELINEERROR_HUMANN} from '../../modules/local/common/pipeline_error'
include { PIPELINEERROR as PIPELINEERROR_KRAKEN2 } from '../../modules/local/common/pipeline_error'
include { PIPELINEEXIT } from '../../modules/local/common/pipeline_exit'

workflow RAPID_TAXONOMIC_PROFILING { 
    take:
    sample_number
    clean_reads      // channel: [ val(id), [ reads1, reads2 ] ]
    samplesheet

    main:
    //check the version of databases
    if(params.metaphlan && params.humann){

        if(params.humann_version =="3"){
            if(!params.mpa_index.contains("vJan21")){
                exit 1, "The version of the MetaPhlAn database compatible with HUMAnN3 is vJan21, please check!"
            }
        }    

        if(params.humann_version == "4"){
            if(!params.mpa_index.contains("vOct22")){
                exit 1, "The version of the MetaPhlAn database compatible with HUMAnN4 is vOct22, please check!"
            }
        }

    }
    
    marker_report = Channel.empty()
    ch_profile = Channel.empty()
    ch_error_log = Channel.empty()

    // run MetaPhlAn or not
    if(params.metaphlan){

        mpa_flag = checkEssentialParams([params.mpa_db])
        if(!mpa_flag) { exit 1, "The required parameter to run MetaPhlAn is: --mpa_db" }
        
        // Create channel for reference databases
        ch_mpa_db = params2Channel(params.mpa_db)

        // metaphlan
        if(params.metaphlan_version == "4.0"){
            METAPHLANV40 (clean_reads, ch_mpa_db)
            ch_mpa_profile = METAPHLANV40.out.profile
            ch_mpa_bowtie2out = METAPHLANV40.out.bowtie2out
        }

        if(params.metaphlan_version == "4.1"){
            METAPHLANV41 (clean_reads, ch_mpa_db)
            ch_mpa_profile = METAPHLANV41.out.profile
            ch_mpa_bowtie2out = METAPHLANV41.out.bowtie2out
        }

        ch_mpa_out = ch_mpa_profile.map{ it-> it[1] }.collect()
        finish_number_mpa = ch_mpa_out.flatten().filter { file -> !file.isEmpty() }.count()

        // merge mpa
        MERGEMPA (sample_number, finish_number_mpa, ch_mpa_out, samplesheet)
        marker_report = MERGEMPA.out.mpa_species
        ch_profile = MERGEMPA.out.mpa_profile
        
        // generate an error log if MetaPhlAn error occurs
        PIPELINEERROR_MPA("MetaPhlAn", sample_number, finish_number_mpa)
        ch_error_log = ch_error_log.mix(PIPELINEERROR_MPA.out.log)

        // sbg2gtdb
        SBG2GTDB (ch_mpa_profile, ch_mpa_db)
        // merge result (to be implemented)
        //MERGESBG2GTDB(SBG2GTDB.out.profile.collect())

        // metaphlan extra abundance compute
        if(params.mpa_extra_abun_method){
            
            ch_mpa_extra_method = Channel.of(params.mpa_extra_abun_method.split(","))
            MPAEXTRAABUN(ch_mpa_extra_method.combine(ch_mpa_bowtie2out), ch_mpa_db)

            MERGEMPAEXTRAABUN(MPAEXTRAABUN.out.profile.groupTuple())   

        }
     
        //run HUMAnN or not
        if(params.humann){
            humann_flag = checkEssentialParams([params.humann_chocophlan_db, params.humann_protein_db, params.humann_map_db])
            if(!humann_flag) { exit 1, "The required parameters to run HUMAnN is:\n --humann_chocophlan_db\n --humann_protein_db\n --humann_map_db" }

            // Create channel for reference databases
            ch_humann_chocophlan_db =  params2Channel(params.humann_chocophlan_db)
            ch_humann_protein_db = params2Channel(params.humann_protein_db)
            ch_humann_map_db =  params2Channel(params.humann_map_db)

            if(params.humann_version == "3"){
                HUMANNV3(clean_reads.join(ch_mpa_profile), ch_humann_chocophlan_db, ch_humann_protein_db, ch_humann_map_db)
                ch_humann_gf = HUMANNV3.out.gf_profile
                ch_humann_ko = HUMANNV3.out.ko_profile
                ch_humann_profile = HUMANNV3.out.profile
            }

            if(params.humann_version == "4"){
              
                taxon4humann = Channel.empty()
                if(params.mpa_options.contains("rel_ab_w_read_stats")){
                    taxon4humann = ch_mpa_profile
                }else{
                    MPAEXTRAABUN4HUMANN(Channel.of("rel_ab_w_read_stats").combine(ch_mpa_bowtie2out), ch_mpa_db)
                    taxon4humann = MPAEXTRAABUN4HUMANN.out.out4humann
                }

                HUMANNV4(clean_reads.join(taxon4humann), ch_humann_chocophlan_db, ch_humann_protein_db, ch_humann_map_db)
                ch_humann_gf = HUMANNV4.out.gf_profile
                ch_humann_ko = HUMANNV4.out.ko_profile
                ch_humann_profile = HUMANNV4.out.profile
            }

            ch_humann_out = ch_humann_gf.map{ it-> it[1] }.collect()
            finish_number_humann = ch_humann_out.flatten().filter { file -> !file.isEmpty() }.count()
            
            MERGEHUMANN(sample_number, finish_number_humann, ch_humann_profile.collect().merge(), samplesheet) 
            ch_profile = ch_profile.mix(MERGEHUMANN.out.humann_profile)
            marker_report = marker_report.mix(MERGEHUMANN.out.humann_report)

            // generate an error log if HUMAnN error occurs
            PIPELINEERROR_HUMANN("HUMAnN", sample_number, finish_number_humann)
            ch_error_log = ch_error_log.mix(PIPELINEERROR_HUMANN.out.log)

            
            def expand_list = [params.humann_map_CARD, params.humann_map_CAZy, params.humann_map_VFDB, params.humann_map_GMMs, params.ko_module, params.ko_pathway].findAll{ it!=null }

            // humann expand CARD CAZY VFDB GMMS
            if(expand_list.size() > 0){
               
                ch_CARD = params.humann_map_CARD ? Channel.fromPath(params.humann_map_CARD).map{ it -> ["CARD", it]} : Channel.empty()
                ch_CAZY = params.humann_map_CAZy ? Channel.fromPath(params.humann_map_CAZy).map{ it -> ["CAZy", it]} : Channel.empty()
                ch_VFDB = params.humann_map_VFDB ? Channel.fromPath(params.humann_map_VFDB).map{ it -> ["VFDB", it]} : Channel.empty()
                ch_GMMS = params.humann_map_GMMs ? Channel.fromPath(params.humann_map_GMMs).map{ it -> ["GMMs", it]} : Channel.empty()               
                ch_gf_expand = ch_CARD.concat(ch_CAZY, ch_VFDB, ch_GMMS)

                ch_module = params.ko_module ? Channel.fromPath(params.ko_module).map{ it -> ["module", it]} : Channel.empty()
                ch_pathway = params.ko_pathway ? Channel.fromPath(params.ko_pathway).map{ it -> ["pathway", it]} : Channel.empty()
                ch_ko_expand = ch_module.concat(ch_pathway)

                // humann expand task number：configured expand database number * samples number
                ch_humann_expand = ch_gf_expand.combine(ch_humann_gf).concat(ch_ko_expand.combine(ch_humann_ko))
                HUMANNEXPAND(ch_humann_expand)

                // group the data by extended database name and combine the profiles of all samples
                MERGEHUMANNEXPAND(HUMANNEXPAND.out.profile.groupTuple(), samplesheet)
                
                //merge all expand profile table 
                ch_expand_profile = MERGEHUMANNEXPAND.out.humann_profile.collectFile(name:"humann.expand.profile.csv")

                ch_profile = ch_profile.mix(ch_expand_profile)

            }
        }
    
    }
    
    // run Kraken2 or not
    if(params.kraken2) {
        
        kraken2_flag = checkEssentialParams([params.kraken2_db])
        if(!kraken2_flag) { exit 1, "The required parameter to run kraken2 is: --kraken2_db" }

        ch_kraken2_db =  params2Channel(params.kraken2_db)

        KRAKEN2(clean_reads, ch_kraken2_db)
        ch_kraken2_out = KRAKEN2.out.species.collect()
        finish_number_kraken2 = ch_kraken2_out.flatten().filter { file -> !file.isEmpty() }.count()

        MERGEKRAKEN2(sample_number, finish_number_kraken2, KRAKEN2.out.mapping.collect(), samplesheet)
        marker_report = marker_report.mix(MERGEKRAKEN2.out.kraken2_species)

        // generate an error log if Kraken2 error occurs
        PIPELINEERROR_KRAKEN2("Kraken2", sample_number, finish_number_kraken2)
        ch_error_log = ch_error_log.mix(PIPELINEERROR_KRAKEN2.out.log)

    }

    // terminate the pipeline if Marker error occurs
    PIPELINEEXIT(ch_error_log.collect())

    marker_report = marker_report.mix(ch_error_log).collect()
    
    emit:
    marker_report
    ch_profile
    
}

