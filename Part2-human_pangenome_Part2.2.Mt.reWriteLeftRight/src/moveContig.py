import os
import re
from Bio import SeqIO
from operator import itemgetter

# 下面的处理方式我暂时不做更改，因为目前来看应该是没有什么问题的，我们可能会有一个 bep的contig对应多个left或者right的序列的情况，但是下面
# 我们写的时候是contig[line[1]] = line[0] 这里key是其他的contig的name，所以就避免了这个问题，同时如果我们存在多匹配的情况，那就自然而然的选择其中的一个。
def pass_contigs(Ensure_contigs, pass_contigs): # 这个是不同的contig的文件。
    F_contigs = open(Ensure_contigs, 'r')
    contigs = {}
    for line in F_contigs.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        if not line == "":
            contigs[line[1]] = line[0]
    F_contigs.close()
    F_contigs = open(pass_contigs, 'r')
    for line in F_contigs.readlines():
        line = line.split('\n')[0]
        line = line.split()
        if not line == "":
            contigs[line[1]] = line[0]
    F_contigs.close()
    return contigs
'''
def cluster_info(cluster_path):  # 这个是representative的序列的位置
    #files = os.listdir(cluster_path)
    contig_id = {}
    for file in files:
        if file[-2:] == "fa":
            for record in SeqIO.parse(cluster_path + file, "fasta"):
                contig_id[str(record.id)] = file.split(".")[0]
    return contig_id
'''

def contig_info(contig_path):  #contig_path: 这个是其他的(maybe is unmapped reads?)
    contigs_seq = {}
    for record in SeqIO.parse(contig_path, "fasta"):
        id = str(record.id).split("=")[0]
        contigs_seq[id] = str(record.seq)
    return contigs_seq





def add_contigs(cluster_path, contigs,contigs_seq):
    f_write = open(cluster_path, 'a')
    for key in contigs:
        #f_write = open(cluster_path + contig_id[contigs[key]] + ".fa", 'a')
        if key in contigs_seq:
            f_write.write(">" + key + "=" +  contigs[key]  + '\n')
            f_write.write(contigs_seq[key]  + '\n')
    f_write.close()



def upRepCluTotUnmap(contigs,othRepTot):
    with open(othRepTot, 'r') as f:
        lines = f.readlines()
    head = 0
    with open(othRepTot, 'w') as f:
        for line in lines:
            if line.strip("\n")[0] == '>':
                header = line[1:(len(line)-1)].split("=")
                if not header[0] in contigs:
                    f.write(line)
                    head = 1
                else:
                    head = 0
            elif head == 1:
                f.write(line)

def upRepCluTotSep(contigs,othRepTot):
    with open(othRepTot, 'r') as f:
        lines = f.readlines()
    head = 0
    with open(othRepTot, 'w') as f:
        for line in lines:
            if line.strip("\n")[0] == '>':
                header = line[1:(len(line)-1)].split("=")
                if not header[0] in contigs and not header[1] in contigs:
                    f.write(line)
                    head = 1
                else:
                    head = 0
            elif head == 1:
                f.write(line)


def moveContig(enSure, Candi, cluster, OtherContig, OtherType):
    # 这里contig key是我们的等待移到其他的文件里面的序列，比如如果是left和unmap来比，那么key就是unmap的序列， value就是能够和这个key比对的left的序列。
    contigs = pass_contigs(enSure, Candi)
    #contig_id = cluster_info(cluster)
    contigs_seq = contig_info(OtherContig)
    # 我们还需要读入
    # 下面的代码在向both的UpRepCluTot里面写的时候就直接写成">加进来的序列=能够和她匹配的bep序列"这种header的形式
    add_contigs(cluster, contigs, contigs_seq)
    root = str(OtherType.split("-")[0])
    othDir = str(OtherType.split("-")[1])
    if re.match(".*unmapped.*",othDir) and len(contigs) != 0 :
        upRepCluTotUnmap(contigs, str(othDir + "/UpG38Unmap.fa"))
    elif len(contigs) != 0 :
        upRepCluTotSep(contigs,str(othDir + "/UpRepCluTot"))
    #这里写的可能有问题，从上下文来看cluster_path应该就是 cluster_folder，就是remain的fasta的地址。



def pass_contigsSepUn(leftDir, rightDir): # 这个是不同的contig的文件。
    leftEn = open(str(leftDir + "/leftUnmapEnsure.txt")) 
    leftCan = open(str(leftDir + "/UpleftUnmapCandi"))

    contigs = {}
    for line in leftEn.readlines():
        line = line.split('\n')[0]
        line = line.split()
        if not line == "":
            #contigs[line[1]] = line[0]
            contigs.setdefault(line[1], [])
            contigs[line[1]].append((line[0],'l'))
    leftEn.close()

    for line in leftCan.readlines():
        line = line.split('\n')[0]
        line = line.split()
        if not line == "":
            #contigs[line[1]] = line[0]
            contigs.setdefault(line[1], [])
            contigs[line[1]].append((line[0],'l'))
    leftCan.close()

    rightEn = open(str(rightDir + "/rightUnmapEnsure.txt"))
    rightCan = open(str(rightDir + "/UprightUnmapCandi"))

    for line in rightEn.readlines():
        line = line.split('\n')[0]
        line = line.split()
        if not line == "":
            print(line)
            #contigs[line[1]] = line[0]
            contigs.setdefault(line[1], [])
            contigs[line[1]].append((line[0],'r'))
    rightEn.close()

    for line in rightCan.readlines():
        line = line.split('\n')[0]
        line = line.split()
        if not line == "":
            #contigs[line[1]] = line[0]
            contigs.setdefault(line[1], [])
            contigs[line[1]].append((line[0],'r'))
    rightCan.close()

    return contigs


def add_contigsSepUn(leftDir, rightDir, contigs,contigs_seq,CompareScore):
    f_writeLeft = open(str(leftDir + "/UpRepCluTot"), 'a')
    f_writeRight = open(str(rightDir + "/UpRepCluTot"), 'a')
    for key in contigs:
        if len(contigs[key]) == 1:
            #f_write = open(cluster_path + contig_id[contigs[key]] + ".fa", 'a')
            if key in contigs_seq:
                if contigs[key][0][1] == "l":
                    f_writeLeft.write(">" + key + "=" +  contigs[key][0][0]  + '\n')
                    f_writeLeft.write(contigs_seq[key]  + '\n')
                else:
                    f_writeRight.write(">" + key + "=" +  contigs[key][0][0]  + '\n')
                    f_writeRight.write(contigs_seq[key]  + '\n')
        else:
            if key in contigs_seq:
                print(key)
                leftAndRight = CompareScore[key]
                leftAndRight = sorted(leftAndRight , key = itemgetter(0,1))
                print(leftAndRight)
                if leftAndRight[1][2] == 'l':
                    f_writeLeft.write(">" + key + "=" +  contigs[key][0][0]  + '\n')
                    f_writeLeft.write(contigs_seq[key]  + '\n')
                else:
                    f_writeRight.write(">" + key + "=" +  contigs[key][1][0]  + '\n')
                    f_writeRight.write(contigs_seq[key]  + '\n')                        

    
def getSepUnCompare(leftDir, rightDir):
    compares = {}
    leftCompare = open(str(leftDir + "/leftUnmapCompare.tsv"),'r')
    rightCompare = open(str(rightDir + "/rightUnmapCompare.tsv"),'r')
    for line in leftCompare.readlines():
        line = line.strip().split()
        compares.setdefault(line[0],[])
        compares[line[0]].append([line[2],(float(line[3])/float(line[4])) + (float(line[3])/float(line[5])), 'l'])
    for line in rightCompare.readlines():
        line = line.strip().split()
        compares.setdefault(line[0],[])
        compares[line[0]].append([line[2],(float(line[3])/float(line[4])) + (float(line[3])/float(line[5])), 'r'])
    leftCompare.close()
    rightCompare.close()
    return compares

def moveContigSepUn(leftDir, rightDir, OtherContig, OtherType):
    # 这里contig key是我们的等待移到其他的文件里面的序列，比如如果是left和unmap来比，那么key就是unmap的序列， value就是能够和这个key比对的left的序列。
    contigs = pass_contigsSepUn(leftDir, rightDir)
    #contig_id = cluster_info(cluster)
    contigs_seq = contig_info(OtherContig)
    CompareScore = getSepUnCompare(leftDir, rightDir)
    # 下面的代码在向both的UpRepCluTot里面写的时候就直接写成">加进来的序列=能够和她匹配的bep序列"这种header的形式
    add_contigsSepUn(leftDir,rightDir, contigs, contigs_seq, CompareScore)

    if re.match(".*unmapped.*",OtherType) and len(contigs) != 0 :
        upRepCluTotUnmap(contigs, str(OtherType + "/UpG38Unmap.fa"))

    #这里写的可能有问题，从上下文来看cluster_path应该就是 cluster_folder，就是remain的fasta的地址。


def deleteSEP(bothDir, SepDir, SepBed):
    contigPos = []
    BothBedFile = open(str(bothDir + "/Both.ending.bed"), 'r')
    SepBedFile = open(SepBed, 'r')
    BothRepUp = open(str(bothDir + "/UpRepCluTot"))
    BothRepUpContigs = []
    SepContigInfo = {}
    removeRep = []
    for lineBothReps in BothRepUp.readlines():
        if lineBothReps[0] == '>':
            lineBothReps = lineBothReps.strip("\n").split("=")
            BothRepUpContigs.append(lineBothReps[0][1:])
    for line in BothBedFile.readlines():
        line = line.strip("\n").split()
        # 这里不用记录both contig的名字，只需要记录位置信息。
        contigPos.append([line[0], line[1], line[2]])
    for lineSep in SepBedFile.readlines():
        lineSep = lineSep.strip("\n").split()
        contigName = str(lineSep[3].split("-")[0] + "-" + lineSep[3].split("-")[1])
        SepContigInfo[contigName] = [lineSep[0], lineSep[1], lineSep[2]]
        # 这里我们或许已经把我们的序列添加到both里面去了。所以我们需要注意我们读取的仍然可能是both的位置
        if contigName in BothRepUpContigs:
            contigPos.append([line[0], line[1], line[2]])

    with open(str(SepDir + "/UpRepCluTot"), 'r') as f:
        lineSepRep = f.readlines()
    headSep = 0

    # 这里有一个包含问题，我们的map文件是总的包含移出去的所有的sep的序列的map信息。而且
    # 这里我们还没有把unmap的序列的信息添加到sep里面，所以肯定是可以都找到对应的位置信息。
    with open(str(SepDir + "/UpRepCluTot"), 'w') as f:
        for lineSepRep in lineSepRep:
            if lineSepRep.strip("\n")[0] == '>':
                contigSep = lineSepRep.strip("\n").split("=")[1] #我这里算的是contig的 represent。
                contigSepInfo = SepContigInfo[contigSep]
                delete = False
                for bothBed in contigPos:
                    if bothBed[0] == contigSepInfo[0]:
                        if( 0 < (int(bothBed[1]) - int(contigSepInfo[2])) <=100)  or ( 0 <(int(contigSepInfo[1]) - int(bothBed[2])) <=100) or ( (int(bothBed[1])<= int(contigSepInfo[1]) <= int(bothBed[2]) ) or (int(bothBed[1]) <= int(contigSepInfo[2]) <= int(bothBed[2]) )):
                            delete = True
                           # print(str(bothBed[0] + "\t" + bothBed[1] + "\t" + bothBed[2] + "\t" + contigSepInfo[0] + "\t" + contigSepInfo[1] + "\t" + contigSepInfo[2] + "\n"))
                if delete == False:
                    f.write(lineSepRep)
                    headSep = 1
                else:
                    headSep = 0
            elif headSep == 1:
                f.write(lineSepRep)


