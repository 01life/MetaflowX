#!/usr/bin/env python

import sys

n = 0
countDir, partentDir, allitem = {},{},[]
with open(sys.argv[1],'r') as inputF, open(sys.argv[2],'w') as outF:
    outF.write("child\tparent\tvalue\n")
    for a in inputF:
        al = a.strip().split('\t')

        if n > 0:
            m = 0
            bl = al[1].strip().split(";")
            for b in bl :
                
                if len(b) > 3:
                    if b not in countDir:
                        countDir[b] = 1
                        allitem.append(b)
                    else:
                        countNuber =  countDir[b] +1
                        countDir[b] = countNuber
                        allitem.append(b)
                else:
                    taxonomayName = b + 'Unclassified'

                if m == 0:
                    partentDir[b] = ','
                else:
                    partentDir[b] = bl[m-1]
                m+=1
        n+=1
    for c in countDir:
            outF.write("%s\t%s\t%s\n"%(c,partentDir[c],str(countDir[c])))
