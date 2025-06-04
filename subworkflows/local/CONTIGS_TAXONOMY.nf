
//
// contigs taxonomy
//

include { params2Channel }              from '../../modules/local/common/utils'
include { checkEssentialParams }        from '../../modules/local/common/utils'
// include { PREBINNING }                  from '../../modules/local/binning/prebinning'
include { KRKEN2CONTIGTAXO }            from '../../modules/local/contig_taxonomy/kraken2_annotate_contigs'
include { ORDERCOMMONID }               from '../../modules/local/contig_taxonomy/order_common_ID'
include { TAXOMETER }                   from '../../modules/local/contig_taxonomy/taxometer'
include { TAXOMETERFILTER }             from '../../modules/local/contig_taxonomy/taxometer_filter'
include { TAXOMETER2BIN }               from '../../modules/local/contig_taxonomy/taxometer2Bin'
include { RENAMECONTIGTAXO }            from '../../modules/local/contig_taxonomy/rename_contig_taxonomy'



workflow CONTIGS_TAXONOMY {

    take:
    contigs           // channel: [ val(id), path(contigs) ]
    depth
    clean_reads       // channel: [ val(id), [ reads1, reads2 ] ]
    contig_taxonomy
    contig_map        // channel: [ val(id), path(contig_map) ]

    main:


        if (params.mode == 5 && params.ContigTaxonomyOptimizer ){

            if (params.contig_kraken2_output && params.contig_taxonomy) {
                contig_taxonomy_input = Channel.fromPath(params.contig_kraken2_output).splitCsv (header:true, sep:',').map { row -> [row.id, row.contig_taxonomy] }
                contig_taxonomy_input.view()
                RENAMECONTIGTAXO(contig_taxonomy_input.join(contig_map))
                ch_contig_taxonomy = RENAMECONTIGTAXO.out.newID_contig_taxonomy


            }else{
                if (!params.contig_taxonomy){
                //step1.1 Kraken2_AnnotateContigs

                    kraken2_flag = checkEssentialParams([params.kraken2_db])
                    if(!kraken2_flag) { exit 1, "The required parameter to run kraken2 is: --kraken2_db" }
                    ch_kraken2_db =  params2Channel(params.kraken2_db)
                    KRKEN2CONTIGTAXO(contigs,ch_kraken2_db)
                    ch_contig_taxonomy =  KRKEN2CONTIGTAXO.out.kraken_taxonomy
                
                }else{
                    exit 1, "[Error] Missing required parameter: --contig_kraken2_output \n\nYou are running MetaflowX in mode 5 using the ContigTaxonomyOptimizer method with `contig_taxonomy = true`. This indicates that the contig taxonomy was obtained in mode 2. Therefore, you must provide the contig taxonomy output file from Kraken2 via the parameter: --contig_kraken2_output.\n\nIf your contigs were not assembled by MetaflowX, you must set `contig_taxonomy = false`. In that case, MetaflowX will automatically classify your contigs.\n\n To resolve this error, either:\n\n\t- Provide the Kraken2 output using `--contig_kraken2_output <file>` <file> Path to a TSV file with 'id' and 'contig_taxonomy' columns, or\n\n\t- Set `contig_taxonomy = false` in your config or command-line arguments."
                }
            }
        }else{
            ch_contig_taxonomy = contig_taxonomy
        }
        


        //step 2 using taxometer improve the taxonomy classification
        ch_contig_taxonomy.view()
        ch_taxometer_input = contigs.join(depth).join(ch_contig_taxonomy)
        ORDERCOMMONID(ch_taxometer_input)
        TAXOMETER(ORDERCOMMONID.out.commoninput)
        TAXOMETERFILTER(TAXOMETER.out.taxometer)

        taximeter95 = TAXOMETERFILTER.out.taxometer95

        TAXOMETER2BIN(contigs.join(TAXOMETERFILTER.out.taxometer95Species))

        //taxometer95Species


        //val(id),path(contigs),path(depth),path(taxonomy)

        //step 3 split to bin

        // refine bin 

        // vamb  recluster  or dastool  comebin



    emit:
    taximeter95

   
}