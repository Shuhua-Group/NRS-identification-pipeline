## iamks 07/26/2022
# It is essential to use this pipeline after pangenome pipeline Part1!
#!/usr/bin/env python

import argparse
import os
import regex as re
from Bio import SeqIO
#import debugpy
import queue
import threading
import logging


# debugpy.listen(6789)
# print("waiting for debug")
# debugpy.wait_for_client()
# print("Attached!")


import sys
import logAndcon
import src.repObtainCluster
import src.moveContig
import src.clusterCompare
import src.mergeLefRight
import src.reorgAlignInof
import src.updataRef
import src.deduplicate
import src.final_ref

def upSEPDatabase(inputDir, IniEnv):
    UpRepTotFile = open(str(inputDir + "/UpRepTot"), 'w')
    for record in SeqIO.parse(str(inputDir + "/UpRepCluTot"), 'fasta'):
        ID = str(record.id).split("=")
        if ID[0] == ID[1]:
            UpRepTotFile.write(">" + ID[0] + "\n")
            UpRepTotFile.write(str(record.seq) + "\n")
    if not os.path.isdir(inputDir + "/UpIndex"): os.mkdir(inputDir + "/UpIndex")
    print(str(IniEnv + " makeblastdb -in " + str(inputDir + "/UpRepTot") + " -dbtype nucl -out " + str(inputDir + "/UpIndex/UpIndex")))
    os.system( IniEnv + " makeblastdb -in " + str(inputDir + "/UpRepTot") + " -dbtype nucl -out " + str(inputDir + "/UpIndex/UpIndex"))




def filter_unmap(novel, unmapDir):
    seq = {}
    UpUnmapFile = open(str(unmapDir + "/UpG38Unmap.fa"), 'w')
    for record in SeqIO.parse(novel, 'fasta'):
        seq[str(record.id)] = [str(record.seq)]
    for record in SeqIO.parse(str(unmapDir + "/Total.unmap"), 'fasta'):
        if record.id in seq:
            UpUnmapFile.write(">" + str(record.id) + "\n")
            UpUnmapFile.write(str(record.seq) + "\n")



def NovelAddSample(filepath, outdir, sample, blastDb):  
    #novelfileOut = open(outPath,"a")
    os.system(str( " module load blast; blastn -query " + filepath  + " -db " + blastDb + " -out " + str(novelDir + "/novel.blast."+ sample) + " -evalue 0.00001 -outfmt 7 -max_target_seqs 1 -num_threads 20")) 
   
    '''
    for line in novelfile.readlines():
        if line[0] == '>':
            line = line.strip()
            novelfileOut.write(str(line + "-" +sample + "\n"))
            sampleID = str(line + "-" +sample)
        else:
            novelfileOut.write(line)
            novel[sampleID] = line.strip()
    return novel


    '''
    # filter the blast result and rewrite the novelSeq
    seqLen = {}
    outFile = open(str(novelDir + "/filter_novel." + sample + ".fa"), 'w')

    for record in SeqIO.parse(filepath, "fasta"):
        seqLen[str(record.id)] = [len(str(record.seq)), str(record.seq)]

    removeContig = []
    with open(str(novelDir + "/novel.blast."+ sample), 'r') as f:
        for line in f.readlines():
            if line[0] != "#":
                line = line.strip().split()
                cov = float(line[3]) / seqLen[line[0]][0]
                if(float(line[2]) > 90 and cov > 0.8):
                    '''
                    outFile.write(str('>' + line[0]) + "\n")
                    outFile.write(seqLen[line[0]][1] + "\n")
                    '''
                    removeContig.append(line[0])
        f.close()

    for seq in seqLen:
        if seq not in removeContig:
            outFile.write(str('>' + seq) + "\n")
            outFile.write(seqLen[seq][1] + "\n")

def multiThreadBlast(q, thread_no, outdir, blastDb):
    while True:
        task = q.get()
        NovelAddSample(task[0], outdir, task[1], blastDb)
        q.task_done()
        print(f'Thread #{thread_no} is doing task #{task} in the queue.')
        
        
     

def RewriteUnmap(UnmapDir, outfile, novel):
    novels ={}
    for record in SeqIO.parse(novel, 'fasta'):
        novels[str(record.id)] = str(record.seq)
    TotalUnmap = open(outfile, 'w')
    for root, dirs, Unmapfiles in os.walk(UnmapDir):
        for files in Unmapfiles:
            if re.match(str("(pan_(.*))-unmap$"), files):
                sample = re.search(str("(pan_(.*))-unmap$"), files).group(1)
                with open(str(root + "/" + files), 'r') as fp:
                    for line in fp:
                        contig = str(line.strip() + "-" + sample)
                        if contig in novels:
                            TotalUnmap.write(">" + contig + "\n")
                            TotalUnmap.write(novels[contig] + "\n")

def filterNovel(novel, Inienv, blastDb,novelDir):
    os.system(str(Inienv + "blastn -query " + novel + " -db " + blastDb + " -out " + str(novelDir + "/novel.blast") + " -evalue 0.00001 -outfmt 7 -max_target_seqs 1 -num_threads 20"))
    seqLen = {}
    outFile = open(str(novelDir + "/filter_novel.fa"), 'w')
    for record in SeqIO.parse(novel, "fasta"):
        seqLen[str(record.id)] = [len(str(record.seq)), str(record.seq)]

    removeContig = []
    with open(str(novelDir + "/novel.blast"), 'r') as f:
        for line in f.readlines():
            if line[0] != "#":
                line = line.strip().split()
                cov = float(line[3]) / seqLen[line[0]][0]
                if(float(line[2]) > 90 and cov > 0.8):
                    '''
                    outFile.write(str('>' + line[0]) + "\n")
                    outFile.write(seqLen[line[0]][1] + "\n")
                    '''
                    removeContig.append(line[0])
        f.close()
    for seq in seqLen:
        if seq not in removeContig:
            outFile.write(str('>' + seq) + "\n")
            outFile.write(seqLen[seq][1] + "\n")

    outFile.close()

def readNovel(filepath,sample):
    novelSeqs = {}
    for record in SeqIO.parse(filepath, 'fasta'):
        novelSeqs[str(record.id) + "-" + sample] = str(record.seq)
    return novelSeqs

def reWriteLeftRight(oriLeft, oriRight, newLeft, newRight):
    leftRes = []
    rightRes = []
    with open(oriLeft, 'r') as fpLeft:
        for line in fpLeft:
            sepLine = line.strip().split()
            if str(sepLine[3]) == "1":
                rightRes.append(str(line))
            else:
                leftRes.append(str(line))
    with open(oriRight, 'r') as fpRight:
        for line in fpRight:
            sepLine = line.strip().split()
            if str(sepLine[3]) == "1":
                leftRes.append(str(line))
            else:
                rightRes.append(str(line))
    with open(newLeft, 'w') as fwLeft:
        for line in leftRes:
            fwLeft.write(line)
    with open(newRight, 'w') as fwRight:
        for line in rightRes:
            fwRight.write(line)

 

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    configPath = "/config"
    parser = argparse.ArgumentParser(description="panGenome pipeline part2")
    parser.add_argument("--indir", help="input dir, pipeline first part output file", required=True, default=None)
    parser.add_argument("--outdir", help="output dir for part2", required=True, default=None)
    parser.add_argument("--config", help="path for configration file ", required=False, default=configPath)
    Args = parser.parse_args()
    if not os.path.isdir(Args.outdir): os.mkdir(Args.outdir)
    config = logAndcon.Config(Args.config)

    loggerFile = str(Args.outdir) + "/pangenomePart2.log"
    logging.basicConfig(filename=str(Args.outdir) + "/pangenomePart2.log2", format="%(asctime)s-%(message)s", level=logging.INFO)
    log = logging.getLogger(__name__)
    
    logger = logAndcon.Logger(loggerFile)

    logger.log(str("inputdir is " + str(Args.indir)))
    logger.log(str("outdir is " + str(Args.outdir)))
    logger.log(str("config is " + str(Args.config)))
    logger.log("READY TO GO!")

    part1Files = {}
    if re.match(".*/$", Args.indir):
        indir = Args.indir
    else:
        indir = Args.indir + "/"
    print(indir)
    leftEndDir = os.path.join(Args.outdir, "leftEnding")
    rightEndDir = os.path.join(Args.outdir, "rightEnding")
    bothEndDir = os.path.join(Args.outdir, "bothEnding")
    unmapDir = os.path.join(Args.outdir, "unmapped")
    novelDir = os.path.join(Args.outdir, "novelSeq")
    if not os.path.isdir(leftEndDir):os.mkdir(leftEndDir)
    if not os.path.isdir(rightEndDir):os.mkdir(rightEndDir)
    if not os.path.isdir(bothEndDir):os.mkdir(bothEndDir)
    if not os.path.isdir(unmapDir):os.mkdir(unmapDir)
    if not os.path.isdir(novelDir): os.mkdir(novelDir)
    IniEnv = str(config.getValue("IniEnv"))
    blastDb  = str(config.getValue("G38path"))

 
    ###########################
    ###########################
    ###########################
    ###########################

    novelFiles = []
    filtedSequence = {}

    
    ###########################
    ###########################
    ###########################
    ###########################

    for root, dirs, files in os.walk(Args.indir, followlinks=True):
        if re.match(str(".*pan_([^/])+$"), root):
            
            sample = root.split("/")[-1]
            logger.log("pre-process for sample : " + sample)

            for file in files:

                #if re.match("unmapped_.*", file):
                #    os.system("cp " + root + "/" + file + " " + unmapDir + "/" + sample + ".unmap")
                #    os.system("cat " + root + "/" + file + " >> " + unmapDir + "/" + "Total.unmap")
                #    part1Files[sample] = [str(unmapDir + "/" + sample + "-unmap")]

                if re.match("novelContig_.*\.fa$", file):
                    filePath = os.path.join(root, file)
                    #novelseq = NovelAddSample(filePath, novelDir, sample)
                    os.system("cp " + root + "/" + file + " " + novelDir + "/" + sample + "-novel")
                    novelFiles.append((filePath,sample))
                    novelseq = readNovel(filePath,sample)
                m = re.match("(.*)MateRegionSort", file)
                if m :
                    os.system("awk -v sample=" + str("pan_" + m.group(1)) + " '{print $0, sample}' " + root + "/" + file + ">>" + Args.outdir + "/mateRegionSort" )

            # os.system("cp " + root + "/contigRegion/leftEndMapping " + leftEndDir + "/" + sample + "-left")
            # os.system("cp " + root + "/contigRegion/rightEndMapping " + rightEndDir + "/" + sample + "-right")

            oriLeft = str(root + "/contigRegion/leftEndMapping")
            oriRight = str(root + "/contigRegion/rightEndMapping")
            newLeft = str(leftEndDir + "/" + sample + "-left")
            newRight = str(rightEndDir + "/" + sample + "-right")

            reWriteLeftRight(oriLeft, oriRight, newLeft, newRight)
            
            logger.log("rewrite left and right over for " + sample)
            
            os.system("cp " + root + "/contigRegion/BothEndMapping " + bothEndDir + "/" + sample + "-both")
            os.system("cp " + root + "/contigRegion/unmappedContig " + unmapDir + "/" + sample + "-unmap")

            logger.log("copy Both and unmap over for " + sample)

            part1Files.setdefault(sample, [])
            part1Files[sample].append(str(leftEndDir + "/" + sample + "-left"))
            part1Files[sample].append(str(rightEndDir + "/" + sample + "-right"))
            part1Files[sample].append(str(bothEndDir + "/" + sample + "-both"))
            #part1Files[sample].append(str(novelDir + "/" + sample + "-novel"))

    forceRewrite = False
    logger.log("copy file over for all samples ! now re write novel sequences")

    if not os.path.exists(str(novelDir + "/filter_novel.fa")) or forceRewrite == True :
        #filterNovel(str(novelDir + "/novel.seq"), IniEnv, blastDb, novelDir)
        q = queue.Queue()

        for i in range(10):
            worker = threading.Thread(target=multiThreadBlast, args=(q, i, novelDir, blastDb), daemon=True)
            worker.start()
        
        for files in novelFiles:
            q.put(files)
        q.join()
        print("blast done for all samples") 
        filterNovel = open(str(novelDir + "/filter_novel.fa"),'w')       
        for root, dirs, novelfiles in os.walk(novelDir):
            for files in novelfiles:
                if re.match(str("filter_novel\.(.*)\.fa"), files):
                    sample = re.search(str("filter_novel\.(.*)\.fa"), files).group(1)
                    
                    for record in SeqIO.parse(str(root + "/" + files),'fasta'):
                        filterNovel.write(">" + str(record.id) + "-" + sample + "\n")
                        filterNovel.write(str(record.seq) + "\n")
                        filtedSequence[(str(record.id) + "-" + sample)] = str(record.seq)
        filterNovel.close()

    logger.log("novel sequence align G38 over, now rewrite unamp sequence")

    RewriteUnmap(unmapDir, str(unmapDir + "/" + "Total.unmap"), str(novelDir + "/filter_novel.fa"))


    print("****************************************get all remapped novel************************************")
    os.system(str(IniEnv + "cd " + bothEndDir + ";awk '{OFS=\"\\t\"} {split(FILENAME,b,\"-\"); if($7==\"1\") print $2,$5,$3,"
                  "$1\"-\"b[1]\"-\"$7\"-\"($8),\"-\";  else print $2,$3,$5,$1\"-\"b[1]\"-\"$7\"-\"($8),\"+\"}' "  + "*-both" " |"
                  " bedtools sort -i > " + Args.outdir + "/bothEnding/Both.ending.bed"))

    os.system(str(IniEnv + "cd " + leftEndDir +";awk '{OFS=\"\\t\"} {split(FILENAME,b,\"-\"); if($4==\"1\") print $2,$7,$8,$1\"-\"b[1]\"-\"$4\"-\"$5\"-\"$3,\"-\"; "
                  " else print $2,$7,$8,$1\"-\"b[1]\"-\"$4\"-\"$5\"-\"$3,\"+\"}' " + "*-left" + " | bedtools sort -i > " + Args.outdir + "/leftEnding/Left.ending.bed"))

    os.system(str(IniEnv + "cd " + rightEndDir + ";awk '{OFS=\"\\t\"} {split(FILENAME,b,\"-\"); if($4==\"1\") print $2,$7,$8,$1\"-\"b[1]\"-\"$4\"-\"$5\"-\"$3,\"-\"; "
                  " else print $2,$7,$8,$1\"-\"b[1]\"-\"$4\"-\"$5\"-\"$3,\"+\"}' " + "*-right" + " | bedtools sort -i > " + Args.outdir + "/rightEnding/Right.ending.bed"))

    os.system(str(IniEnv + "bedtools merge -d 20 -c 4 -o distinct -i " + Args.outdir + "/bothEnding/Both.ending.bed > " + Args.outdir + "/bothEnding/merge_both.bed "))
    os.system(str(IniEnv + "bedtools merge -d 20 -c 4 -o distinct -i " + Args.outdir + "/leftEnding/Left.ending.bed > " + Args.outdir + "/leftEnding/merge_left.bed "))
    os.system(str(IniEnv + "bedtools merge -d 20 -c 4 -o distinct -i " + Args.outdir + "/rightEnding/Right.ending.bed > " + Args.outdir + "/rightEnding/merge_right.bed "))


    
    ###########################
    ###########################
    ###########################
    ###########################
    
    if not os.path.isdir(leftEndDir + "/leftRep"):os.mkdir(leftEndDir + "/leftRep")
    if not os.path.isdir(rightEndDir + "/rightRep"):os.mkdir(rightEndDir + "/rightRep")
    if not os.path.isdir(bothEndDir + "/bothRep"):os.mkdir(bothEndDir + "/bothRep")
    if not os.path.isdir(leftEndDir + "/leftCluster"):os.mkdir(leftEndDir + "/leftCluster")
    if not os.path.isdir(rightEndDir + "/rightCluster"):os.mkdir(rightEndDir + "/rightCluster")
    if not os.path.isdir(bothEndDir + "/bothCluster"):os.mkdir(bothEndDir + "/bothCluster")


    leftRep, leftCluster  = str(leftEndDir + "/leftRep/"), str(leftEndDir + "/leftCluster/")
    rightRep, rightCluster = str(rightEndDir + "/rightRep/"), str(rightEndDir + "/rightCluster/")
    bothRep, bothCluster = str(bothEndDir + "/bothRep/"), str(bothEndDir + "/bothCluster/")


    ###########################
    ###########################
    ###########################
    ###########################


    
    src.repObtainCluster.rep_cluster(str(novelDir + "/filter_novel.fa"), str(leftEndDir + "/merge_left.bed"), leftRep, leftCluster, str(leftEndDir + "/leftFullCluster"))
    src.repObtainCluster.rep_cluster(str(novelDir + "/filter_novel.fa"), str(rightEndDir + "/merge_right.bed"), rightRep, rightCluster,str(rightEndDir + "/rightFullCluster"))
    src.repObtainCluster.rep_cluster(str(novelDir + "/filter_novel.fa"), str(bothEndDir + "/merge_both.bed"), bothRep, bothCluster, str(bothEndDir + "/bothFullCluster"))
    filter_unmap(str(novelDir + "/filter_novel.fa"), unmapDir)


    src.repObtainCluster.nucmerRep(str(leftRep), str(leftCluster),IniEnv)
    src.repObtainCluster.nucmerRep(str(rightRep), str(rightCluster),IniEnv)
    src.repObtainCluster.nucmerRep(str(bothRep), str(bothCluster),IniEnv)

    if not os.path.isdir(leftEndDir + "/leftIndex"): os.mkdir(leftEndDir + "/leftIndex")
    if not os.path.isdir(rightEndDir + "/rightIndex"): os.mkdir(rightEndDir + "/rightIndex")
    if not os.path.isdir(bothEndDir + "/bothIndex"): os.mkdir(bothEndDir + "/bothIndex")
    os.system(IniEnv + " makeblastdb -in " + leftEndDir + "/RepTotal -dbtype nucl -out " + str(leftEndDir + "/leftIndex/leftindex"))
    os.system(IniEnv + " makeblastdb -in " + rightEndDir + "/RepTotal -dbtype nucl -out " + str(rightEndDir + "/rightIndex/rightindex"))
    os.system(IniEnv + " makeblastdb -in " + bothEndDir + "/RepTotal -dbtype nucl -out " + str(bothEndDir + "/bothIndex/bothindex"))


    # compair lep / rep / unplaced with Bep
    os.system(IniEnv + "blastn -db " + str(bothEndDir + "/bothIndex/bothindex ") + "-query " + str(leftEndDir + "/leftFullCluster ") + " "
        "-outfmt \"6  qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapopen gaps evalue bitscore\" -max_target_seqs 1 "
        " -max_hsps 1 -out " +  str(bothEndDir + "/bothLeftCompare.tsv"))

    os.system(IniEnv + "blastn -db " + str(bothEndDir + "/bothIndex/bothindex ") + "-query " + str(rightEndDir + "/rightFullCluster ") + " "
        "-outfmt \"6  qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapopen gaps evalue bitscore\" -max_target_seqs 1 "
        " -max_hsps 1 -out " +  str(bothEndDir + "/bothRightCompare.tsv"))

    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$4>=0.99 && ($6-$13) /$5>=0.9 ) print $2,$1}' " + str(bothEndDir + "/bothLeftCompare.tsv") + " > " + str(bothEndDir + "/bothLeftEnsure.txt"))
    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$5<0.9 && ($6-$13)/$4>=0.99 && ($6-$13)/$5>=0.8) print $2,$1}' " + str(bothEndDir + "/bothLeftCompare.tsv") + " > " + str(bothEndDir + "/bothLeftCandidate.txt"))

    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$4>=0.99 && ($6-$13) /$5>=0.9 ) print $2,$1}' " + str(bothEndDir + "/bothRightCompare.tsv") + " > " + str(bothEndDir + "/bothRightEnsure.txt"))
    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$5<0.9 && ($6-$13)/$4>=0.99 && ($6-$13)/$5>=0.8) print $2,$1}' " + str(bothEndDir + "/bothRightCompare.tsv") + " > " + str(bothEndDir + "/bothRightCandidate.txt"))


    os.system(IniEnv + "blastn -db " + str(bothEndDir + "/bothIndex/bothindex ") + "-query " + str(unmapDir +  "/" + "UpG38Unmap.fa") + " "
        "-outfmt \"6  qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapopen gaps evalue bitscore\" -max_target_seqs 1 "
        " -max_hsps 1 -out " +  str(bothEndDir + "/bothUnmapCompare.tsv"))

    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$4>=0.99 && ($6-$13) /$5>=0.9 ) print $2,$1}' " + str(bothEndDir + "/bothUnmapCompare.tsv") + " > " + str(bothEndDir + "/bothUnmapEnsure.txt"))
    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$5<0.9 && ($6-$13)/$4>=0.99 && ($6-$13)/$5>=0.85 ) print $2,$1}' " + str(bothEndDir + "/bothUnmapCompare.tsv") + " > " + str(bothEndDir + "/bothUnmapCandidate.txt"))


    src.clusterCompare.EnsureCandidate(str(bothEndDir + "/bothLeftCandidate.txt"), Args.outdir, str(bothEndDir + "/Both.ending.bed"), str(bothEndDir + "/UpBothLeftCandi"))
    src.clusterCompare.EnsureCandidate(str(bothEndDir + "/bothRightCandidate.txt"), Args.outdir, str(bothEndDir + "/Both.ending.bed"), str(bothEndDir + "/UpBothRightCandi"))
    src.clusterCompare.EnsureCandidate(str(bothEndDir + "/bothUnmapCandidate.txt"), Args.outdir, str(bothEndDir + "/Both.ending.bed"), str(bothEndDir + "/UpBothUnmapCandi"))


    src.moveContig.moveContig(str(bothEndDir + "/bothLeftEnsure.txt"), str(bothEndDir + "/UpBothLeftCandi"),str(bothEndDir + "/UpRepCluTot"), str(leftEndDir + "/UpRepCluTot"),str(bothEndDir + "-" + leftEndDir))
    src.moveContig.moveContig(str(bothEndDir + "/bothRightEnsure.txt"), str(bothEndDir + "/UpBothRightCandi"),str(bothEndDir + "/UpRepCluTot"), str(rightEndDir + "/UpRepCluTot"),str(bothEndDir + "-" + rightEndDir))

    print("before we delete all single end sequence within 100bp of bothEnding sequence ")
    src.moveContig.deleteSEP(bothEndDir, leftEndDir, str(leftEndDir + "/Left.ending.bed"))
    src.moveContig.deleteSEP(bothEndDir, rightEndDir, str(rightEndDir + "/Right.ending.bed"))
    
    src.moveContig.moveContig(str(bothEndDir + "/bothUnmapEnsure.txt"), str(bothEndDir + "/UpBothUnmapCandi"),str(bothEndDir + "/UpRepCluTot"), str(unmapDir + "/UpG38Unmap.fa"),str(bothEndDir + "-" + unmapDir))


    # up data our database !
    upSEPDatabase(leftEndDir,IniEnv)
    upSEPDatabase(rightEndDir, IniEnv)
    print("move for both is over --------------- exit")
    os.system(IniEnv + "blastn -db " + str(leftEndDir + "/UpIndex/UpIndex ") + "-query " + str(unmapDir +  "/" + "UpG38Unmap.fa") + " "
        "-outfmt \"6  qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapopen gaps evalue bitscore\" -max_target_seqs 1 "
        " -max_hsps 1 -out " +  str(leftEndDir + "/leftUnmapCompare.tsv"))

    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$4>=0.99 && ($6-$13) /$5>=0.9 ) print $2,$1}' " + str(leftEndDir + "/leftUnmapCompare.tsv") + " > " + str(leftEndDir + "/leftUnmapEnsure.txt"))
    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$5<0.9 && ($6-$13)/$4>=0.99 && ($6-$13)/$5>=0.85) print $2,$1}' " + str(leftEndDir + "/leftUnmapCompare.tsv") + " > " + str(leftEndDir + "/leftUnmapCandidate.txt"))


    src.clusterCompare.EnsureCandidate(str(leftEndDir + "/leftUnmapCandidate.txt"), Args.outdir,str(leftEndDir + "/Left.ending.bed"), str(leftEndDir + "/UpleftUnmapCandi"))
    #src.moveContig.moveContig(str(leftEndDir + "/leftUnmapEnsure.txt"), str(leftEndDir + "/UpleftUnmapCandi"),str(leftEndDir + "/UpRepCluTot"), str(unmapDir + "/UpG38Unmap.fa"),str(leftEndDir + "-" + unmapDir))


    os.system(IniEnv + "blastn -db " + str(rightEndDir + "/UpIndex/UpIndex ") + "-query " + str(unmapDir +  "/" + "UpG38Unmap.fa") + " "
        "-outfmt \"6  qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapopen gaps evalue bitscore\" -max_target_seqs 1 "
        " -max_hsps 1 -out " +  str(rightEndDir + "/rightUnmapCompare.tsv"))

    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$4>=0.99 && ($6-$13) /$5>=0.9 ) print $2,$1}' " + str(rightEndDir + "/rightUnmapCompare.tsv") + " > " + str(rightEndDir + "/rightUnmapEnsure.txt"))
    os.system("awk '{OFS=\"\\t\"}{if($3>99 && ($6-$13)/$5<0.9 && ($6-$13)/$4>=0.99 && ($6-$13)/$5>=0.85 ) print $2,$1}' " + str(rightEndDir + "/rightUnmapCompare.tsv") + " > " + str(rightEndDir + "/rightUnmapCandidate.txt"))


    src.clusterCompare.EnsureCandidate(str(rightEndDir + "/rightUnmapCandidate.txt"), Args.outdir,str(rightEndDir + "/Right.ending.bed"), str(rightEndDir + "/UprightUnmapCandi"))
    print("***********************unmapped read had compared with left and right **************************")
    
  
    ###########################
    ###########################
    ###########################
    ###########################

    src.moveContig.moveContigSepUn(leftEndDir , rightEndDir , str(unmapDir + "/UpG38Unmap.fa"),unmapDir)
    

    

    if not os.path.isdir(Args.outdir + "/MergeLefRig"): os.mkdir(Args.outdir + "/MergeLefRig")


    src.mergeLefRight.Merge(leftEndDir, rightEndDir, str(Args.outdir + "/MergeLefRig"), config)

    #paths = [str(Args.outdir + "/MergeLefRig/identity"), str(Args.outdir + "/MergeLefRig/contain"),str(Args.outdir + "/MergeLefRig/overlap"),str(Args.outdir + "/MergeLefRig/part")]
   



    src.reorgAlignInof.alignment_priority(str(leftEndDir+ "/Left.ending.bed"), str(rightEndDir+ "/Right.ending.bed"), str(Args.outdir + "/MergeLefRig"))


    src.reorgAlignInof.redefine_alignment_result(str(Args.outdir + "/MergeLefRig/identity"), str(Args.outdir + "/MergeLefRig/"))
    if not os.path.isdir(Args.outdir + "/MergeLefRig/PartConfirm"):os.mkdir(Args.outdir + "/MergeLefRig/PartConfirm")


    # confirm the part align
    src.reorgAlignInof.conFirmPart(str(Args.outdir + "/MergeLefRig"), config)
    


    # merge and updata
    src.reorgAlignInof.PopinsMerge(str(Args.outdir + "/MergeLefRig"), config)
    


    src.updataRef.upDataRef(Args.outdir)


    '''
    ###########################
    ###########################
    ###########################

    src.updataRef.mergeAndAlign(Args.outdir, config)
        



    src.deduplicate.rmRedundance(str(Args.outdir + "/AllPlaced/AllPlacedAligned.tsv"), 2000, 98, 0.95, Args.outdir)
 

    # remove redundancy for unplaced contig
    src.final_ref.finalRef(Args.outdir)




    '''
    ###########################
    ###########################
    ###########################
