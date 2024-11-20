#!/usr/bin/env python

import argparse

def markerQS(completeness: float, contamination: float) -> str:
    """
    Determine quality score level based on completeness and contamination.
    Reference: https://www.nature.com/articles/nbt.3893/tables/1
    """
    if completeness > 90 and contamination < 5:
        return 'HQ'
    elif completeness >= 50 and contamination < 10:
        return 'MQ'
    else:
        return 'LQ'

def parse_normal_QS(original_QSFile: str) -> dict:
    tmpDir = {}
    with open(original_QSFile, 'r') as original_QSF:
        header = original_QSF.readline()
        for line in original_QSF:
            al = line.strip().split("\t")
            qslevel = markerQS(float(al[2]), float(al[3]))
            tmpDir[al[0]] = [al[2], al[3], qslevel]
    return tmpDir

def parse_ReAss_QS(reass_QSFile: str) -> dict:
    tmpDir = {}
    with open(reass_QSFile, 'r') as reass_QSF:
        header = reass_QSF.readline()
        for line in reass_QSF:
            bl = line.strip().split("\t")
            binid = bl[0].split("_reassembly_contigs")[0]
            qslevel = markerQS(float(bl[1]), float(bl[2]))
            tmpDir[binid] = [bl[1], bl[2], qslevel]
    return tmpDir

def parse_ReBin_QS(reBin_QSFile: str) -> dict:
    tmpDir = {}
    with open(reBin_QSFile, 'r') as rebin_QSF:
        header = rebin_QSF.readline()
        for line in rebin_QSF:
            cl = line.strip().split("\t")
            qslevel = markerQS(float(cl[1]), float(cl[2]))
            tmpDir[cl[0]] = [cl[1], cl[2], qslevel]
    return tmpDir

def parse_ReRefine_QS(reRefine_QSFile: str) -> dict:
    tmpDir = {}
    with open(reRefine_QSFile, 'r') as refine_QSF:
        for line in refine_QSF:
            dl = line.strip().split("\t")
            qslevel = markerQS(float(dl[2]), float(dl[3]))
            tmpDir[dl[0]] = [dl[2], dl[3], qslevel]
    return tmpDir

def parse_improve_info(improveFile: str) -> dict:
    tmpDir = {}
    with open(improveFile, 'r') as improve_QSF:
        for line in improve_QSF:
            el = line.strip().split("\t")
            qslevel = markerQS(float(el[1]), float(el[2]))
            tmpDir[el[0]] = [el[1], el[2], qslevel, el[3], el[4]]
    return tmpDir

def get_level_and_quality(data_dict, onebin, default_level, default_quality):
    if onebin in data_dict:
        level = data_dict[onebin][-1]
        quality = f"{onebin}\t{data_dict[onebin][0]}\t{data_dict[onebin][1]}"
    else:
        level = default_level
        quality = ''
    return level, quality




def main(pickbin_file, original_QS_file,reAss_QS_file, reBin_QS_file, reRefine_QS_file,improveFile, output_file,plot_prefix):
    original_qs_dir = parse_normal_QS(original_QS_file)
    reass_qs_dir 	= parse_ReAss_QS(reAss_QS_file)
    rebin_qs_dir 	= parse_ReBin_QS(reBin_QS_file)
    refine_qs_dir 	= parse_ReRefine_QS(reRefine_QS_file)
    improve_dir 	= parse_improve_info(improveFile)


    pick_list = []

    with open(pickbin_file,'r') as pickbin_F, open(output_file, 'w') as output_F, open(plot_prefix + ".sanket.txt",'w') as sankey_F, open(plot_prefix + ".point.txt",'w') as point_F,open(plot_prefix + "_final_improvement.xls",'w') as improvement_F:
        improvement_F.write(f"BinID\tCompleteness\tContamination\tLevel\tStaus\tStage\n")
        sankey_F.write(f"BinID\tInit\tReAss\tReBin\tReFine\tFinal\n")
        point_F.write(f"BinID\tCompleteness\tContamination\tStage\n")
        output_F.write(f"binid\tinit_Completeness\tinit_Contamination\tinit_QS_Level\tafter_reass_Completeness\tafter_reass_Contamination\tafter_reass_QS_Level\tafter_rebin_Completeness\tafter_rebin_Contamination\tiafter_rebin_QS_Level\tafter_refine_Completeness\tafter_refine_Contamination\tafter_refine_QS_Level\n")
        for  p in pickbin_F:
            pl  = p.strip().split("\t")
            pickbin = pl[0]
            pick_list.append(pickbin)
            outtxt = '\t'.join(original_qs_dir[pickbin])
            outtxt += '\t' + '\t'.join(reass_qs_dir.get(pickbin, ['-', '-', '-']))
            outtxt += '\t' + '\t'.join(rebin_qs_dir.get(pickbin, ['-', '-', '-']))
            outtxt += '\t' + '\t'.join(refine_qs_dir.get(pickbin, ['-', '-', '-']))
            output_F.write(f"{pickbin}\t{outtxt}\n")

            improvement_F.write(str(pickbin)+"\t"+'\t'.join(improve_dir.get(pickbin, original_qs_dir[pickbin] + ['failed', 'Init']))+"\n")

        for onebin in pick_list:
            if onebin not in original_qs_dir:
                continue

            init_level, init_quality = get_level_and_quality(original_qs_dir, onebin, '', '')
            ass_level, ass_quality = get_level_and_quality(reass_qs_dir, onebin, init_level, init_quality)
            bin_level, bin_quality = get_level_and_quality(rebin_qs_dir, onebin, ass_level, ass_quality)
            refine_level, refine_quality = get_level_and_quality(refine_qs_dir, onebin, bin_level, bin_quality)

            sankey_F.write(f"{onebin}\t{init_level}\t{ass_level}\t{bin_level}\t{refine_level}\t{refine_level}\n")

            for stage, quality in [('Init', init_quality), ('ReAss', ass_quality), 
                                ('ReBin', bin_quality), ('ReFine', refine_quality)]:
                if quality != "":
                    point_F.write(f"{quality}\t{stage}\n")
                else:
                    print(f"remove{onebin}\t{stage}")
                    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some bin files.")
    parser.add_argument('--pickbin', type=str, required=True, help='Path to the pick bin information txt')
    parser.add_argument('--original_QS', type=str, required=True, help='Path to the original bin QS file')
    parser.add_argument('--ReAss_QS', type=str, required=True, help='Path to the reassembly QS file')
    parser.add_argument('--ReBin_QS', type=str, required=True, help='Path to the rebin QS file')
    parser.add_argument('--ReRefine_QS', type=str, required=True, help='Path to the rerefine QS file')
    parser.add_argument('--improve', type=str, required=True, help='Path to the improve info file')
    parser.add_argument('--outfile', type=str, required=True, help='Output file')
    parser.add_argument('--platout', type=str, required=True, help='plat output prefix')
    args = parser.parse_args()
    main(args.pickbin,args.original_QS, args.ReAss_QS, args.ReBin_QS, args.ReRefine_QS, args.improve, args.outfile, args.platout)