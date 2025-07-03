
//
// contigs taxonomy
//

include { params2Channel }                                              from '../../modules/local/common/utils'
include { checkEssentialParams }                                        from '../../modules/local/common/utils'
include { KRKEN2CONTIGTAXO }                                            from '../../modules/local/contig_taxonomy/kraken2_annotate_contigs'
include { CATCONTIG }                                                   from '../../modules/local/contig_taxonomy/cat_contig'

include { ORDERCOMMONID }                                               from '../../modules/local/contig_taxonomy/order_common_ID'
// include { TAXOMETER }                                                   from '../../modules/local/contig_taxonomy/taxometer'


include { TAXVAMB }                                                     from '../../modules/local/contig_taxonomy/taxvamb'
include { TAXVAMB2BIN }                                                 from '../../modules/local/contig_taxonomy/taxvamb2Bin'
include { RENAMECONTIGTAXO }                                            from '../../modules/local/contig_taxonomy/rename_contig_taxonomy'

include { TAXVAMBFILTER }                                               from '../../modules/local/contig_taxonomy/taxvamb_filter'

include { TAXOMETERFILTER }                                             from '../../modules/local/contig_taxonomy/taxometer_filter'
include { TAXOMETER2BIN }                                               from '../../modules/local/contig_taxonomy/taxometer2Bin'


workflow CONTIGS_TAXONOMY {

    take:
    contigs           // channel: [ val(id), path(contigs) ]
    depth
    clean_reads       // channel: [ val(id), [ reads1, reads2 ] ]
    contig_taxonomy
    contig_map        // channel: [ val(id), path(contig_map) ]
    prodigal_faa      // channel: [ val(id), path(faa) ]

    main:


        if (params.mode == 5 && params.ContigTaxonomyOptimizer ){

            if (params.contig_taxonomy_output && params.contig_taxonomy) {
                contig_taxonomy_input = Channel.fromPath(params.contig_taxonomy_output).splitCsv (header:true, sep:',').map { row -> [row.id, row.contig_taxonomy] }
                contig_taxonomy_input.view()
                RENAMECONTIGTAXO(contig_taxonomy_input.join(contig_map))
                ch_contig_taxonomy = RENAMECONTIGTAXO.out.newID_contig_taxonomy


            }else{
                if (!params.contig_taxonomy){


                    if (params.contig_taxonomy && !(params.kraken2_contig || params.cat_contig)) {exit 1, "When enabling contig taxonomy analysis (--contig_taxonomy), you must set either --kraken2_contig or --cat_contig to true."}

                    //step1.1 Kraken2_AnnotateContigs
                     if (params.kraken2_contig){
                        kraken2_flag = checkEssentialParams([params.kraken2_db])
                        if(!kraken2_flag) { exit 1, "The required parameter to run kraken2 is: --kraken2_db" }
                        ch_kraken2_db =  params2Channel(params.kraken2_db)
                        KRKEN2CONTIGTAXO(contigs,ch_kraken2_db)
                        ch_contig_taxonomy =  KRKEN2CONTIGTAXO.out.kraken_taxonomy
                    }

                     if (params.cat_contig){
                        ch_cat_pack = params2Channel(params.cat_pack)
                        ch_cat_db = params2Channel(params.cat_gtdb_db)
                        CATCONTIG(contigs.join(prodigal_faa),ch_cat_db,ch_cat_pack)
                        ch_contig_taxonomy = CATCONTIG.out.cat_taxonomy
                    }

                
                }else{
                    exit 1, "[Error] Missing required parameter: --contig_taxonomy_output \n\nYou are running MetaflowX in mode 5 using the ContigTaxonomyOptimizer method with `contig_taxonomy = true`. This indicates that the contig taxonomy was obtained in mode 2. Therefore, you must provide the contig taxonomy output file from Kraken2 or CAT via the parameter: --contig_taxonomy_output.\n\nIf your contigs were not assembled by MetaflowX, you must set `contig_taxonomy = false`. In that case, MetaflowX will automatically classify your contigs.\n\n To resolve this error, either:\n\n\t- Provide the Kraken2 or CAT output using `--contig_taxonomy_output <file>` <file> Path to a TSV file with 'id' and 'contig_taxonomy' columns, or\n\n\t- Set `contig_taxonomy = false` in your config or command-line arguments."
                }
            }
        }else{
            ch_contig_taxonomy = contig_taxonomy
        }
        


        //step 2 using taxometer improve the taxonomy classification
        ch_contig_taxonomy.view()
        ch_taxometer_input = contigs.join(depth).join(ch_contig_taxonomy)
        
        //TaxVamb
        ORDERCOMMONID(ch_taxometer_input)
        TAXVAMB(ORDERCOMMONID.out.commoninput)


        // if not using checkm to filter the bin quality, need to change the code of TAXVAMB.
        // ch_checkm2_db = params2Channel(params.checkm2_db)
        // TAXVAMBFILTER(TAXVAMB.out.bins, ch_checkm2_db)



        //Taxometer
        
        TAXOMETERFILTER(TAXVAMB.out.taxometer)
        taximeter95 = TAXOMETERFILTER.out.taxometer95
        taxometer95deepestLevel = TAXOMETERFILTER.out.taxometer95deepest


        TAXOMETER2BIN(contigs.join(taxometer95deepestLevel))

        

        //val(id),path(contigs),path(depth),path(taxonomy)

        //step 3 split to bin

        // refine bin 

        // vamb  recluster  or dastool  comebin



    emit:
    taximeter95

   
}