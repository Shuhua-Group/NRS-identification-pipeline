import os
import collections
from Bio import SeqIO
import regex as re
import sys
import logging
import argparse
#import debugpy

# debugpy.listen(6789)
# print("wait for debug")
# debugpy.wait_for_client()
# print("attached")
 
def insertion_points(LEP_folder, REP_folder):

    contig_insertion = {}
    files = os.listdir(LEP_folder)
    for file in files:
        if re.match(str(".*-left$"), file):
            f = open(LEP_folder + "/" + file, 'r')
            m = re.match("(.*)-left", file)
            sample_ID = m[1]
            for line in f.readlines():
                line = line.split('\n')[0]
                line = line.split()
                contig_ID = line[0] + '-' + sample_ID

                length = int(line[4])
                #if line[3] == "1":
                # 我清楚为什么要用101了，因为之前这里它选择的map的contig 末端的长度是100！！！！！所以如果是reverse，他要用101
                #    map_len = 150 - int(line[8])
                #else:
                map_len = int(line[8]) 
                # 我来解释一下，下面 line[7]是我们的ref的end，也就是末端比对在ref上的end位置。
                # 可能是因为在它输出的时候没有输出contig的长度，所以它建立了fai index 然后通过读取index来或者contig的
                # 这里有有点没理解，这里第一个应该是insertion point， 但是左右的insertion point是不一样的，这里的选取是不是有点问题。
                contig_insertion[contig_ID] = [int(line[7]) - 1, 1.0 *  int(line[8]) / length, length, map_len]
            f.close()
    files = os.listdir(REP_folder)
    for file in files:
        if re.match(str(".*-right$"), file):
            f = open(REP_folder + "/" + file, 'r')
            m = re.match("(.*)-right", file)
            sample_ID = m[1]
            for line in f.readlines():
                line = line.split('\n')[0]
                line = line.split('\t')
                contig_ID = line[0] + '-' + sample_ID

                length = int(line[4])
               # if line[3] == "1":
                # 我清楚为什么要用101了，因为之前这里它选择的map的contig 末端的长度是100！！！！！所以如果是reverse，他要用101
               #     map_len = 151 - int(line[8])
               # else:
                map_len = int(line[8])
                # 我来解释一下，下面 line[7]是我们的ref的end，也就是末端比对在ref上的end位置。
                # 可能是因为在它输出的时候没有输出contig的长度，所以它建立了fai index 然后通过读取index来或者contig的
                # 
                contig_insertion[contig_ID] = [int(line[7]) - 1, 1.0 *  int(line[8]) / length, length,
                                           map_len]
            f.close()

    return contig_insertion


def obtain_seq(UpLeftTot, UpRightTot, UpBothTot):
    LEP_rep = collections.defaultdict(list)
    LEP_cluster = collections.defaultdict(list)
    REP_rep = collections.defaultdict(list)
    REP_cluster = collections.defaultdict(list)
    Both_cluster = collections.defaultdict(list)
    # cluster
   # files = os.listdir(LEP_cluster_folder)
   # for file in files:
   # 这个name其实就是contig name，因为之前我们写cluster和representative的contig的时候我们是分序列来写的所以contig有
   #     name = file.split(".")[0]
    for record in SeqIO.parse(UpLeftTot, "fasta"):
        leftCluContig = str(record.id).split("=")[0]
        leftRepContig = str(record.id).split("=")[1]
        if leftCluContig == leftRepContig:
            LEP_rep[leftRepContig].append(str(record.seq))
            LEP_cluster[leftRepContig].append([record.id, str(record.seq)])
        else:
            LEP_cluster[leftRepContig].append([record.id, str(record.seq)])

    for record in SeqIO.parse(UpRightTot, "fasta"):
        rightCluContig = str(record.id).split("=")[0]
        rightRepContig = str(record.id).split("=")[1]
        if rightCluContig == rightRepContig:
            REP_rep[rightRepContig].append(str(record.seq))
            REP_cluster[rightRepContig].append([record.id, str(record.seq)])
        else:
            REP_cluster[rightRepContig].append([record.id, str(record.seq)])

    for record in SeqIO.parse(UpBothTot, "fasta"):
        bothCluContig = str(record.id).split("=")[0]
        bothRepContig = str(record.id).split("=")[1]
        if bothCluContig == bothRepContig:
            #REP_rep[rightRepContig].append(str(record.seq))
            Both_cluster[bothRepContig].append([record.id, str(record.seq)])
        else:
            Both_cluster[bothRepContig].append([record.id, str(record.seq)])


    return LEP_cluster, REP_cluster, Both_cluster , LEP_rep, REP_rep

# 后面4个中前俩是colloction 后面俩是set ，后面四个都是空的，第一个文件是我们的left和right比对的结果，也就是 part；identity 等等的文件。
# 然后LEP和REP 我们的contig的list，也就是文件名对应序列head和sequence的dictionary。
# file是我们输入的文件，然后LEP_cluster以及REP_cluster是我们的lep以及rep序列的cluster的文件读入header以及到一个list里面的结果
# 这里我们做的事情就是查找，在clust的文件里面查找


# 这个函数的作用就是
def move_SEP_contigs(file, LEP_cluster, REP_cluster, LEP_add, REP_add, delete_LEP, delete_REP, BEP_add):
    filename = file.split("/")[-1]
    f = open(file, 'r')

    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split()

        rep = line[len(line) - 1].split('-')
        # 我们当时写的时候是在 left right left.l 或者 如果是right是representative的话就是left right right.r
        if rep[len(rep) - 1] == 'l':
            #下面是把left和right进行了合并。
            for number1 in range(len(LEP_cluster[line[0]])):
                firstID = LEP_cluster[line[0]][number1][0].split("=")[0]
                if filename == "confirm.final_part.txt":
                    LEP_add[str(rep[0] + "-" + rep[1])].append([str(firstID + "=" + rep[0] + "-" + rep[1]), LEP_cluster[line[0]][number1][1]])
                else:
                    BEP_add[str(rep[0] + "-" + rep[1])].append([str(firstID + "=" + rep[0] + "-" + rep[1]), LEP_cluster[line[0]][number1][1]])

            for number1 in range(len(REP_cluster[line[1]])):
                firstID = REP_cluster[line[1]][number1][0].split("=")[0]
                if filename == "confirm.final_part.txt":
                    LEP_add[str(rep[0] + "-" + rep[1])].append([str(firstID + "=" + rep[0] + "-" + rep[1]), REP_cluster[line[1]][number1][1]])
                else:
                    BEP_add[str(rep[0] + "-" + rep[1])].append([str(firstID + "=" + rep[0] + "-" + rep[1]), REP_cluster[line[1]][number1][1]])

        elif rep[len(rep) - 1] == 'r':
            for number1 in range(len(LEP_cluster[line[0]])):
                firstID = LEP_cluster[line[0]][number1][0].split("=")[0]
                if filename == "confirm.final_part.txt":
                    REP_add[str(rep[0] + "-" + rep[1])].append([str(firstID + "=" + rep[0] + "-" + rep[1]), LEP_cluster[line[0]][number1][1]])
                else:
                    BEP_add[str(rep[0] + "-" + rep[1])].append([str(firstID + "=" + rep[0] + "-" + rep[1]), LEP_cluster[line[0]][number1][1]])
            for number1 in range(len(REP_cluster[line[1]])):
                firstID = REP_cluster[line[1]][number1][0].split("=")[0]
                if filename == "confirm.final_part.txt":
                    REP_add[str(rep[0] + "-" + rep[1])].append([str(firstID + "=" + rep[0] + "-" + rep[1]), REP_cluster[line[1]][number1][1]])
                else:
                    BEP_add[str(rep[0] + "-" + rep[1])].append([str(firstID + "=" + rep[0] + "-" + rep[1]), REP_cluster[line[1]][number1][1]])

        delete_LEP.add(line[0])
        delete_REP.add(line[1])
    print("left delete identiy size is " + str(len(delete_LEP)))
    print("right delete identiy size is " + str(len(delete_REP)))
    f.close()

    return LEP_add, REP_add, delete_LEP, delete_REP , BEP_add



def is_file_all_whitespace(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if not line.isspace():
                return False
    return True

# 这个函数要返回一个seq的名称和序列以及对应的left， right insertion point， 对应的left和right的contig map site。
# 也就是经过检查，合格的序列会放到这个list里面。
# info = [seq1, seq2, LeftInsertion, rightInsertion, leftPoint, rightPoint, f"mergeOver-{contig_ID}", sequence, repSeq[1], "LsmR"(或者RsmL)]
# 上面我就没过滤left inserte point 小于 right insert point 的情况。
# transBed[contig_ID] = info
# 这个函数可以用来找overlap ， 也可以用来找contain， 都是ok的。如果是contain的话，就用最后的那个repSeq[1]来判断这条序列是比对到left还是right。
def overlap_SEP(file, contig_insert_pos, ArgOut, mergeDir):
    f = open(file, 'r')
    move_LEP = set()
    move_REP = set()
    transBed ={}
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split()
        seq1 = line[0]  # LEP
        seq2 = line[1]  # REP
        contig_ID = seq1 + "-" + seq2
        repSeq = line[2].split("-")
        sample1 = seq1.split("-")[1]
        contig1 = seq1.split("-")[0]
        sample2 = seq2.split("-")[1]
        contig2 = seq2.split("-")[0]
        with open(str(ArgOut + "/leftEnding/" + sample1 + "-left"), 'r') as fp:
            for line in fp:
                line = line.strip().split()
                if line[0] == contig1:
                    LeftcontigMapPos = int(line[8])
                    LeftcontigLen = int(line[4])
                    LeftReverse = line[3]
                    LeftInsertion = line[2]
                    LeftChr = line[1]

        with open(str(ArgOut + "/rightEnding/" + sample2 + "-right"), 'r') as fp:
            for line in fp:
                line = line.strip().split()
                if line[0] == contig2:
                    RightcontigMapPos = int(line[8])
                    RightcontigLen = int(line[4])
                    RightReverse = line[3]
                    RightInsertion = line[2]
                    RightChr = line[1]
        mergedSeqDir = mergeDir + "/popinsMerge/" + contig_ID
        mergedSeq = mergedSeqDir + "/" + contig_ID + ".fa"
        if os.path.exists(mergedSeq) and not is_file_all_whitespace(mergedSeq):
            seq1Fa = mergedSeqDir + "/" + seq1 + "/contigs.fa"
            seq2Fa = mergedSeqDir + "/" + seq2 + "/contigs.fa"
            cmd1 = "module load mummer; nucmer -p " + mergedSeqDir + "/left.align " + seq1Fa + " " + mergedSeq
            cmd2 = "module load mummer; nucmer -p " + mergedSeqDir + "/right.align " + seq2Fa + " " + mergedSeq
            print(cmd1 + "\n")
            print(cmd2 + "\n")
            os1 = os.system(cmd1)
            os2 = os.system(cmd2)
            if os1 ==0 and os2 == 0:
                print("align merge successful")
            else:
                print("align merge failed, program stopped")
                sys.exit(1)
            leftAlign = mergedSeqDir + "/left.align.delta"
            rightAlign = mergedSeqDir + "/right.align.delta"
            leftAlignRefStart = None
            leftAlignRefEnd = None
            leftAlignConStart = None
            leftAlignConEnd = None
            rightAlignRefStart = None
            rightAlignRefEnd = None
            rightAlignConStart = None
            rightAlignConEnd = None
            leftPoint = None 
            rightPoint = None
            with open(leftAlign , 'r') as fp:
                for line in fp:
                    if line[0] != ">":
                        line = line.strip().split()
                        if len(line) == 7:
                            leftAlignConStart = int(line[0])
                            leftAlignConEnd = int(line[1])
                            leftAlignRefStart = int(line[2])
                            leftAlignRefEnd = int(line[3])
                            if LeftReverse == "0":
                                leftPoint = LeftcontigMapPos - leftAlignConStart + leftAlignRefStart
                            else:
                                newPoint = LeftcontigLen - (LeftcontigMapPos - 1)
                                leftPoint = newPoint - leftAlignConStart + leftAlignRefStart
            if(leftPoint == None):
                print("get merge point for leftside of  overlap contig failed, abort!")
                sys.exit()  
            else:
                print(f"leftPoint for mergeseq {contig_ID} is " + str(leftPoint))
            

            with open(rightAlign , 'r') as fp:
                for line in fp:
                    if line[0] != ">":
                        line = line.strip().split()
                        if len(line) == 7:
                            rightAlignConStart = int(line[0])
                            rightAlignConEnd = int(line[1])
                            rightAlignRefStart = int(line[2])
                            rightAlignRefEnd = int(line[3])
                            if RightReverse == "0":
                                rightPoint = RightcontigMapPos - rightAlignConStart + rightAlignRefStart
                            else:
                                newPoint = RightcontigLen - (RightcontigMapPos - 1)
                                rightPoint = newPoint - rightAlignConStart + rightAlignRefStart
            if(rightPoint == None):
                print("get merge point for rightside of  overlap contig failed, abort!")
                sys.exit()  
            else:
                print(f"leftPoint for mergeseq {contig_ID} is " + str(rightPoint))

            for record in SeqIO.parse(mergedSeq, 'fasta'):
                sequence = str(record.seq)
            
            if leftPoint < rightPoint:
                #transBed.setdefault(f"{seq1}-{seq2}", [])
                info = [seq1, seq2, LeftInsertion, RightInsertion, leftPoint, rightPoint, f"{contig_ID}-merge", sequence, repSeq[2], "LsmR", LeftChr, RightChr]
                transBed[contig_ID] = info
            else:
                info = [seq1, seq2, LeftInsertion, RightInsertion, leftPoint, rightPoint, f"{contig_ID}-merge", sequence, repSeq[2], "RsmL", LeftChr, RightChr]
                transBed[contig_ID] = info
        else:
            print(f"merge failed for {seq1} and {seq2}, do not move these two sequences.")



    f.close()
    return transBed





# transBed 里面的信息。 info = [seqName1, seqName2, LeftInsertion, rightInsertion, leftPoint, rightPoint, f"{contig_ID}-merge", MergedSequence, repSeq[1]]
def check_BEP_contigs(file , transBed, LEP_cluster, REP_cluster, delete_LEP, delete_REP, LEP_add, REP_add,
                      LEP_rep, REP_rep, merged_contig_folder,BEP_add, contain):
    f = open(file, 'r')
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split()
        LEP_ID = line[0]
        REP_ID = line[1]
        rep = LEP_ID + "-" + REP_ID
        repSeq = line[len(line)-1]
        if repSeq[len(repSeq) - 1] == "l":
            repSEQ = line[0]
        else:
            repSEQ = line[1]
        #for record in SeqIO.parse(merged_contig_folder +  "/" + rep + "/" + rep + '.fa', "fasta"):
        #    sequence = str(record.seq)
        if rep in transBed:
            info = transBed[rep]
            if contain == True:
                if repSeq[len(repSeq) - 1] == "l":
                    # LEP_rep 这个变量基本没啥用。
                    LEP_rep[rep] = [info[7]]
                    LEP_add[str(info[6])].append([str(info[6]), info[7]])
                    for number1 in range(len(LEP_cluster[line[0]])):
                        # 不写 cluster中的rep
                        seqName =  LEP_cluster[line[0]][number1][0].split("=")[0] 
                        if not LEP_cluster[line[0]][number1][0].split("=")[0] == LEP_cluster[line[0]][number1][0].split("=")[1]:
                            LEP_add[str(info[6])].append([f"{seqName}={info[6]}", LEP_cluster[line[0]][number1][1]])
                    for number1 in range(len(REP_cluster[line[1]])):
                        seqName =  REP_cluster[line[1]][number1][0].split("=")[0] 
                        if not REP_cluster[line[1]][number1][0].split("=")[0] == REP_cluster[line[1]][number1][0].split("=")[1]:
                            LEP_add[str(info[6])].append([f"{seqName}={info[6]}", REP_cluster[line[1]][number1][1]])
                elif repSeq[len(repSeq) - 1] == "r":
                    REP_rep[rep] = [info[7]]
                    REP_add[str(info[6])].append([str(info[6]), info[7]])
                    # 这里如果按照本文件来看，我们的一个cluster的文件对应的是一个list，list里面是header和sequence
                    for number1 in range(len(LEP_cluster[line[0]])):
                        seqName = LEP_cluster[line[0]][number1][0].split("=")[0] 
                        if not LEP_cluster[line[0]][number1][0].split("=")[0] == LEP_cluster[line[0]][number1][0].split("=")[1]:
                            REP_add[str(info[6])].append([f"{seqName}={info[6]}", LEP_cluster[line[0]][number1][1]])
                    for number1 in range(len(REP_cluster[line[1]])):
                        seqName = REP_cluster[line[1]][number1][0].split("=")[0] 
                        if not REP_cluster[line[1]][number1][0].split("=")[0] == REP_cluster[line[1]][number1][0].split("=")[1]:
                            REP_add[str(info[6])].append([f"{seqName}={info[6]}", REP_cluster[line[1]][number1][1]])
                delete_LEP.add(line[0])
                delete_REP.add(line[1])
            else:
                if info[9] == "LsmR":
                    BEP_add[info[6]].append([str(info[6]), info[7]])
                    for number1 in range(len(LEP_cluster[line[0]])):
                        seqName = LEP_cluster[line[0]][number1][0].split("=")[0] 
                        if not LEP_cluster[line[0]][number1][0].split("=")[0] == LEP_cluster[line[0]][number1][0].split("=")[1]:
                            BEP_add[info[6]].append([f"{seqName}={info[6]}", LEP_cluster[line[0]][number1][1]])
                    for number1 in range(len(REP_cluster[line[1]])):
                        seqName = REP_cluster[line[1]][number1][0].split("=")[0] 
                        if not REP_cluster[line[1]][number1][0].split("=")[0] == REP_cluster[line[1]][number1][0].split("=")[1]:
                            BEP_add[info[6]].append([f"{seqName}={info[6]}", REP_cluster[line[1]][number1][1]])

                    delete_LEP.add(line[0])
                    delete_REP.add(line[1])

    f.close()
    return BEP_add, LEP_rep, REP_rep, LEP_add, REP_add, delete_LEP, delete_REP


def reobtain_rep_cluster(LEP_cluster, REP_cluster, LEP_rep, REP_rep, BEP_add, LEP_add, REP_add, delete_LEP, delete_REP,Both_cluster,
                         rightUpdata, leftUpdata, bothUpdata, merged_contig_folder, transBedOver, transBedCon):
    right_rep_path = open(rightUpdata + '.fa', 'w')
    testSumRight = 0
    for key in REP_cluster:
        if key not in delete_REP:
            for number in range(len(REP_cluster[key])):
                testSumRight = testSumRight +1 
                right_rep_path.write(">" + REP_cluster[key][number][0] + "\n")
                right_rep_path.write(REP_cluster[key][number][1] + "\n")
    print("the writen rep sum for right is " + str(testSumRight))

    left_rep_path = open(leftUpdata + '.fa', 'w')
    testSumLeft = 0
    for key in LEP_cluster:
        if key not in delete_LEP:
            for number in range(len(LEP_cluster[key])):
                testSumLeft = testSumLeft +1 
                left_rep_path.write(">" + LEP_cluster[key][number][0] + "\n")
                left_rep_path.write(LEP_cluster[key][number][1] + "\n")
    print("the writen rep sum for left is " + str(testSumLeft))

    for key in REP_add:
        for number in range(len(REP_add[key])):
            right_rep_path.write(">" + REP_add[key][number][0] + "\n")
            right_rep_path.write(REP_add[key][number][1] + "\n")



    for key in LEP_add:
        for number in range(len(LEP_add[key])):
            left_rep_path.write(">" + LEP_add[key][number][0] + "\n")
            left_rep_path.write(LEP_add[key][number][1] + "\n")

    left_rep_path.close()
    right_rep_path.close()
    with open(str(bothUpdata + '.fa'), 'w') as both_rep_path:
        # add 2EP contigs
        for key in Both_cluster:
            for number in range(len(Both_cluster[key])):
                 both_rep_path.write(str(">" + Both_cluster[key][number][0]) + "\n")
                 both_rep_path.write(str(Both_cluster[key][number][1] ) + "\n")

        for key in BEP_add:
            '''
            for record in SeqIO.parse(merged_contig_folder + key + "/" + key + '.fa', "fasta"):
                both_rep_path.write('>' + key + "=" + key + '\n')
                both_rep_path.write(str(record.seq) + '\n')
            '''
            for number in range(len(BEP_add[key])):
                both_rep_path.write(str(">" + BEP_add[key][number][0] + "\n"))
                both_rep_path.write(str(BEP_add[key][number][1] + "\n"))
    both_rep_path.close()
    # transBed 里面的信息。 info = [seqName1, seqName2, LeftInsertion, rightInsertion, leftPoint, rightPoint, f"{contig_ID}-merge", MergedSequence, repSeq[1], "LsmR", LeftChr, RightChr]
    with open(merged_contig_folder + "/MergeOverlap.position", 'w') as fw:
        fw.write("ContigName\tChr\tLeftInsertPoint\tLeftContigPoint\tRightInsertPoint\tRightContigPoint\tReverse\tRepSeq\n")
        for contigID in transBedOver:
            info = transBedOver[contigID]
            fw.write(f"{info[6]}\t{info[10]}\t{info[2]}\t{info[4]}\t{info[3]}\t{info[5]}\t0\t{info[8]}\n")

    with open(merged_contig_folder + "/MergeContain.position", 'w') as fw:
        fw.write("ContigName\tChr\tLeftInsertPoint\tLeftContigPoint\tRightInsertPoint\tRightContigPoint\tReverse\tRepSeq\n")
        for contigID in transBedCon:
            info = transBedCon[contigID]
            fw.write(f"{info[6]}\t{info[10]}\t{info[2]}\t{info[4]}\t{info[3]}\t{info[5]}\t0\t{info[8]}\n")

    
def upDataRef(ArgOut):

    # 在对代码有一定的了解以后这里也不难理解，首先就是我们的left 或者 right的存储形式是在文件夹里面有以rep序列名字命名的文件。
    # 然后下面我们读取的时候是按照rep序列的名称和它对应的rep序列或者cluster序列来读取的
    LEP_cluster, REP_cluster, Both_cluster, LEP_rep, REP_rep = obtain_seq(str(ArgOut + "/leftEnding/UpRepCluTot"), str(ArgOut + "/rightEnding/UpRepCluTot"), str(ArgOut + "/bothEnding/UpRepCluTot"))

    # 上面这个函数就是把所有的contig 以及rep的序列都读到list里面，这里首先我们得清楚我们需要的文件是有两部分。left的和right的。
    # 具体的就是用差不多类似一个字典，字典的key是文件名，value是序列的头和sequence所组成的list。
    alignment_files = ["final_identity.txt", "confirm.final_part.txt"]
    LEP_add = collections.defaultdict(list)
    REP_add = collections.defaultdict(list)
    BEP_add = collections.defaultdict(list)
    delete_LEP = set()
    delete_REP = set()
    for file in alignment_files:
        file_path = ArgOut + "/MergeLefRig/" + file
        LEP_add, REP_add, delete_LEP, delete_REP, BEP_add = move_SEP_contigs(file_path, LEP_cluster, REP_cluster, LEP_add,
                                                                    REP_add, delete_LEP, delete_REP, BEP_add)
    file = ArgOut + "/MergeLefRig/" + "final_overlap.txt"

    # 我个人的感觉，这个insertion point的作用是把insert的情况有给写了一下，这里的contig fai明显是rep或者lep里面
    # 的contig的index，因为我们从这个里面提取信息，另外这里的rep和lep应该就是没有经过awk转化成为bed的原始的文件
    # insertion_points这个函数的作用差不多就是把我们的rep或者lep里面的序列的insertion信息给提取然后拓展一下。

    contig_insertion = insertion_points(str(ArgOut + "/leftEnding"), str(ArgOut + "/rightEnding"))
    # 然后我们根据上一行判断的insert的信息然后判断我们的overlap的那个文件的结果，通过结合从而进行筛选。
    # 我们的move做的是什么呢？我们的move会对overlap里面的序列进行判断，其中一部份不符合要求的序列分类到--移动到SEP的类别，其他的才能移动到BEP，算是一个粗分类。
    transBedOver = overlap_SEP(file, contig_insertion, str(ArgOut), str(ArgOut)+"/MergeLefRig")

    # check_BEP_contigs函数算是接着上面的来，对粗分类的结果进一步的划分到不同的add里面，方便我们最终书写，add就是rep以及对应的cluster序列，cluster序列从cluster的collection中提取
    # 下面这个函数是怎么做的呢？？下面这个函数的做法是我们从overlap的这些情况里面判断BEP是否存在。我们判断的标准就是这个contig是不是在
    # move_REP或者move_LEP之类的东西里面存在，如果不存在呢么可能是一个BEP
    # 其实整个过程是一个筛选的过程，把overlap里面的东西按照不同的标准划分到 lep rep以及bep里面，也就是 BEP_add, LEP_add, REP_add,
    # 同时也更新了 REP/LEP rep的东西，也是根据overlap里面的内容来更新的，然后所有的东西，不论判断标准怎么都写到delete里面
    # 当然上面更新 rep的时候用的是merge的东西来更新的
    
    fileCon = ArgOut + "/MergeLefRig/final_contained.txt"
    transBedCon = overlap_SEP(fileCon, contig_insertion, str(ArgOut), str(ArgOut)+"/MergeLefRig")

    BEP_add, LEP_rep, REP_rep, LEP_add, REP_add, delete_LEP, delete_REP = check_BEP_contigs(file, transBedOver,
                                                                                            LEP_cluster, REP_cluster,
                                                                                            delete_LEP, delete_REP,
                                                                                            LEP_add, REP_add, LEP_rep,
                                                                                            REP_rep,str(ArgOut + "/MergeLefRig/popinsMerge"), BEP_add, False)
    BEP_add, LEP_rep, REP_rep, LEP_add, REP_add, delete_LEP, delete_REP = check_BEP_contigs(fileCon, transBedCon,
                                                                                            LEP_cluster, REP_cluster,
                                                                                            delete_LEP, delete_REP,
                                                                                            LEP_add, REP_add, LEP_rep,
                                                                                            REP_rep,str(ArgOut + "/MergeLefRig/popinsMerge"),BEP_add, True)
    reobtain_rep_cluster(LEP_cluster, REP_cluster, LEP_rep, REP_rep, BEP_add, LEP_add, REP_add, delete_LEP, delete_REP,Both_cluster ,
                         str(ArgOut + "/rightEnding/UpRepCluTot2"), str(ArgOut + "/leftEnding/UpRepCluTot2"),  str(ArgOut + "/bothEnding/UpRepCluTot2") , str(ArgOut + "/MergeLefRig"), transBedOver, transBedCon)


def mergeAndAlign(ArgsOutdir, config):
    IniEnv = str(config.getValue("IniEnv"))
    if not os.path.isdir(ArgsOutdir + "/AllPlaced"):os.mkdir(ArgsOutdir + "/AllPlaced")

    AllPlaced = open(str(ArgsOutdir + "/AllPlaced/AllplacedContig.fa"),'w')
    with open(str(ArgsOutdir + "/rightEnding/UpRepCluTot2.fa"), 'r') as f:
        for line in f.readlines():
            AllPlaced.write(line)
    f.close()
    with open(str(ArgsOutdir + "/leftEnding/UpRepCluTot2.fa"), 'r') as f:
        for line in f.readlines():
            AllPlaced.write(line)
    f.close()
    with open(str(ArgsOutdir + "/bothEnding/UpRepCluTot2.fa"), 'r') as f:
        for line in f.readlines():
            AllPlaced.write(line)
    f.close()
    AllPlaced.close()
    if not os.path.isdir(ArgsOutdir + "/AllPlaced/dbindex"): os.mkdir(ArgsOutdir + "/AllPlaced/dbindex")
    os.system(str(IniEnv + "makeblastdb -in " + str(ArgsOutdir + "/AllPlaced/AllplacedContig.fa") + " -dbtype nucl -out " + ArgsOutdir + "/AllPlaced/dbindex/dbindex"))
    os.system(str(IniEnv + "blastn -db " + ArgsOutdir + "/AllPlaced/dbindex/dbindex -query " + ArgsOutdir + "/AllPlaced/AllplacedContig.fa" + " -outfmt \"6  qseqid sseqid  pident slen qlen length qstart qend sstart send mismatch gapopen gaps evalue bitscore \"  -out "
                  + ArgsOutdir + "/AllPlaced/AllPlacedAligned.tsv"))