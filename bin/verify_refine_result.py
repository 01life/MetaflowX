import sys

orgQS, refineQS = {},{}
refinebinPath = {}

with open(sys.argv[1],'r') as orignal_QSF, open(sys.argv[2],'r') as refineQSF, open(sys.argv[3],'r')  as bin_pathF, open(sys.argv[4],'r') as mapF, open(sys.argv[5],'w') as outF:

    for a in orignal_QSF:
        al = a.strip().split("\t")
        orgQS[al[0]] = al

    for b in refineQSF:
        bl = b.strip().split("\t")
        binid = bl[0].lstrip("Refine_")
        refineQS[binid] = bl

    for c in bin_pathF:
        cl  =  c.strip().split("\t")
        refinebinPath[cl[0]] = cl[1]

    
    header =  mapF.readline()
    outF.write(header)

    for d in mapF:
        binID,binGenome,refGenomes,highestAbundanceFq1,fq1 = d.strip().split("\t")
        if binID in refineQS  and binID in orgQS:
            orgCompleteness = orgQS[binID][2]
            orgContamination = orgQS[binID][3]
            refineCompleteness = refineQS[binID][1]
            refineContamination = refineQS[binID][2]

            if float(refineCompleteness) >= float(orgCompleteness) or float(refineContamination) <= float(orgContamination):
                binGenome = refinebinPath[binID]
                
            outF.write(f"{binID}\t{binGenome}\t{refGenomes}\t{highestAbundanceFq1}\t{fq1}\n")





    




