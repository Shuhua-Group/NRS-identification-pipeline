from pickle import TRUE
from tkinter import FALSE
import regex as re
import os
from Bio import SeqIO

def Merge(leftDir, rightDir, outDir, config):

    IniEnv = str(config.getValue("IniEnv"))
    leftPos = {}
    rightPos = {}
    matchPair = {}
    leftRep = []
    rightRep = []

    #four situation of result:
    identity = []
    overlap = []
    contain = []
    part = []


    with open(str(leftDir + "/UpRepCluTot"), 'r') as f:
        for line in f:
            if line[0] == '>':
                line = line.strip("\n").split("=")
                leftRep.append(line[1])

    with open(str(rightDir + "/UpRepCluTot"), 'r') as f:
        for line in f:
            if line[0] == '>':
                line = line.strip("\n").split("=")
                rightRep.append(line[1])


    with open(str(leftDir + "/Left.ending.bed"),'r') as f:
        for line in f:
            line = line.strip("\n").split()
            contigName = line[3].split("-")[0] + "-" + line[3].split("-")[1]
            Insertion = line[3].split("-")[4]
            if contigName in leftRep:
                leftPos[contigName] = [line[0], line[1],line[2], line[4], Insertion]

    with open(str(rightDir + "/Right.ending.bed"),'r') as f:
        for line in f:
            line = line.strip("\n").split()
            contigName = line[3].split("-")[0] + "-" + line[3].split("-")[1]
            Insertion = line[3].split("-")[4]
            if contigName in rightRep:
                rightPos[contigName] = [line[0], line[1], line[2], line[4], Insertion]
                for leftContig in leftPos:
                    if leftPos[leftContig][0] == line[0] and max(abs(int(line[1]) - int(leftPos[leftContig][1])), abs(int(line[1]) - int(leftPos[leftContig][2])), abs(int(line[2]) - int (leftPos[leftContig][1])), abs(int(line[2]) - int (leftPos[leftContig][2]))) <= 100 :
                        matchPair.setdefault(leftContig, [])
                        matchPair[leftContig].append(contigName)


    leftSeq ={}

    for record in SeqIO.parse(str(leftDir + "/UpRepCluTot"), "fasta"):
        clusterRep = record.id.strip("\n").split("=")[1]
        contigname = record.id.strip("\n").split("=")[0]
        leftSeq[(contigname, clusterRep)] = str(record.seq)

    rightSeq = {}

    for record in SeqIO.parse(str(rightDir + "/UpRepCluTot"), "fasta"):
        clusterRep = record.id.strip("\n").split("=")[1]
        contigname = record.id.strip("\n").split("=")[0]
        rightSeq[(contigname, clusterRep)] = str(record.seq)

    for left in matchPair:
        #right = matchPair[left]
        for right in matchPair[left]:
            outDirMatch = str(outDir + "/" + left + "-" + right)
            if not os.path.isdir(outDirMatch): os.mkdir(outDirMatch)
            leftRep = outDirMatch + "/" + left + ".fa"
            rightRep = outDirMatch + "/" + right + ".fa"
            leftReverse = leftPos[left][3]
            rightReverse = rightPos[right][3]
            leftInsert = leftPos[left][4]
            rightInsert = rightPos[right][4]


            with open(leftRep, 'w') as f:
                f.write(">" + left + "\n")
                f.write(leftSeq[left, left] + "\n")

            with open(rightRep, 'w') as f:
                f.write(">" + right + "\n")
                f.write(rightSeq[right, right] + "\n")

            leftCluster = outDirMatch + "/" + left + ".cluster.fa"
            rightCluster = outDirMatch + "/" + right + ".cluster.fa"

            with open(leftCluster, 'w') as f:
                for i in leftSeq:
                    if i[1] == left:
                        f.write(">" + i[0] + "\n")
                        f.write(leftSeq[i] + "\n")

            with open(rightCluster, 'w') as f:
                for i in rightSeq:
                    if i[1] == right:
                        f.write(">" + i[0] + "\n")
                        f.write(rightSeq[i] + "\n")
        
            os.system(IniEnv + "nucmer -p " + outDirMatch + "/alignRes " + leftRep + " " + rightRep)
            os.system(IniEnv + "delta-filter -q -r -g -m -1 " + outDirMatch + "/alignRes.delta > " + outDirMatch + "/filter_Rep_delta" )
            os.system(IniEnv + "show-coords -H -T -l -c -o " + outDirMatch + "/filter_Rep_delta > " + outDirMatch + "/filter_Rep_coords" )
            with open(str(outDirMatch + "/filter_Rep_coords"), 'r') as f:
                print("get filted coords" + str(outDirMatch + "/filter_Rep_coords") )
                #classed = False
                for line in f:
                    line = line.strip("\n").split()
                    if float(line[6]) >= 97 and float(line[9]) >=70 and  float(line[10])>=70: #and classed == False:
                        add = TRUE
                        for exiden in identity:
                            if (str(exiden[1][11] + exiden[1][12]) == str(line[11] + line[12])):
                                if (float(exiden[1][6]) <= float(line[6])):
                                    identity.remove(exiden)
                                    identity.append((outDirMatch, line))
                                add = FALSE
                        if add == TRUE:
                            identity.append((outDirMatch, line))       
                    elif (float(line[6]) >= 90) and float(line[9])>=60 and float(line[10])>=60: #and classed == False :
                        if (float(leftInsert) <= float(rightInsert)):
                            overlap.append((outDirMatch, line))
                        else:
                            contain.append((outDirMatch, line ))
                    elif (float(line[6]) >= 50) and float(line[9]) >= 60 and float(line[10])>=60: #and classed == False:
                        part.append((outDirMatch, line))

    IdentityFile = str(outDir + "/identity")
    OverlapFile = str(outDir + "/overlap")
    ContainFile = str(outDir + "/contain")
    PartFile = str(outDir + "/part")

    with open(IdentityFile,'w') as f:
        for i in identity:
            line = '\t'.join(i[1])
            f.write(line + "\n")
    with open(OverlapFile,'w') as f:
        for i in overlap:
            line = '\t'.join(i[1])
            f.write(line + "\n")
    with open(ContainFile,'w') as f:
        for i in contain:
            line = '\t'.join(i[1])
            f.write(line + "\n")
    with open(PartFile,'w') as f:
        for i in part:
            line = '\t'.join(i[1])
            f.write(line + "\n")

