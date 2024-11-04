#!/usr/bin/env python

import sys,numpy

parent_child_Dir , childValueDir, allrelationship, allitem, passitem= {},{},{},[],[]

with open(sys.argv[1],'r') as mpaF, open(sys.argv[2],'w') as outF:
    outF.write("child\tparent\tgrandsonNumber\tabundance\n")
    for a in mpaF:
        al = a.strip().split('\t')
        idline = al[0].split("|")
        if len(idline) > 1 and a[0] in ["d","k"]:
                parent_child_Dir.setdefault(idline[-2],[]).append(idline[-1])
                avgAbundance = numpy.average(list(map(float, al[1:])))
                childValueDir[idline[-1]] = [idline[-2],avgAbundance] # Sum the total abundance of all samples.
                allitem.append(idline[-1])
                allrelationship[idline[-1]] = idline

        else:
            if a[0] in ["d","k"] :
                avgAbundance = numpy.average(list(map(float, al[1:])))
                #childValueDir[al[0]] = ["",al[-1]]
                childValueDir[al[0]] = ["",avgAbundance]
                allitem.append(al[0])


    uniqitem = set(allitem)
    
    #check each item should has parent and parent has value.
    for c in uniqitem:
        n = 0
        if c[0] not in ["d","k"] :
            for d in allrelationship[c]:
                if d not in uniqitem:
                    n += 1
                if n > 0:
                    if d not in passitem:
                        passitem.append(d)
    print(passitem)
    for b in uniqitem:
        outputline = ''
        if b[0] in ["d","k"] :
            outputline = "%s\t%s\t%s\t%s\n"%(b,",",str(len(parent_child_Dir[b])),str(childValueDir[b][1]))
        else:
            parent = childValueDir[b][0]
            m = 0
            for e in allrelationship[b]:
                if e in passitem:
                    m += 1
            if m < 1:
                if b in parent_child_Dir:
                    outputline = "%s\t%s\t%s\t%s\n"%(b,childValueDir[b][0],str(len(parent_child_Dir[b])),str(childValueDir[b][1]))
                else:
                    if b[0] == "s":
                        outputline = "%s\t%s\t%s\t%s\n"%(b,childValueDir[b][0],"1",str(childValueDir[b][1]))
            else:
                print(str(b)+"\t: there is some family member was abandoned \tpass")
        if  outputline != "":
            outF.write(str(outputline))
        else:
            print(str(b)+"\t: did not get info")