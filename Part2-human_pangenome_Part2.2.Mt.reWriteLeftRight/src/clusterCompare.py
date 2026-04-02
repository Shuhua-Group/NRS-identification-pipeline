import re
import logging
log = logging.getLogger(__name__)

#src.clusterCompare.EnsureCandidate(str(bothEndDir + "/bothLeftCandidate.txt"), Args.outdir, str(bothEndDir + "/Both.ending.bed"), str(bothEndDir + "/UpBothLeftCandi"))
def EnsureCandidate(CandiFile, rootPath, bedFile, outfile): # bed file both或者single end的bed文件，有其在染色体上的位置信息。
    BedFile = open(bedFile, 'r')
    bedInfo = {}
    for lineBed in BedFile.readlines():
        lineBed = lineBed.strip().split()
        contigINFO = lineBed[3].split("-")
        bedInfo[str(contigINFO[0] + "-" + contigINFO[1])] = [lineBed[0], lineBed[1], lineBed[2]]
    candiFile = open(CandiFile, 'r')
    outFile = open(outfile, 'a')
    MateLinkFile = open(str(rootPath + "/mateRegionSort"), 'r')
    MateLinkDirc = {}
    with open(str(rootPath + "/mateRegionSort"), 'r') as fp:
        for line in fp:
            line = line.strip().split()
            contigName = line[1] + "-" + line[len(line)-1]
            log.info(f"anchor-info-found-for-{contigName}")
            MateLinkDirc.setdefault(contigName, [])
            MateLinkDirc[contigName].append([line[6],line[8],line[9]])

    for line in candiFile.readlines():
        
        line = line.strip().split()
        queryContigName = line[1]
        refContigName = line[0]
        #contigName = contigName.split('-')
        #m = re.match("pan_(.*)", contigName[1])
        #sampleName = m.gourp(1)
        #contig = contigName[0]
        log.info(f"ReadCandi-processing-for-query-{queryContigName}-and-ref-{refContigName}")
        if refContigName not in bedInfo:
            continue
        # 这里要思考的一个问题就是有没有，我们用bed文件来 直接取某一个可以和unmapp比对上的left序列的位置，能不能肯定的得到结果
        # 这里要分两种情况来考虑，第一种是我们的BEP文件和left，right，unmap比对后的确认，第二种是 我们的SEP和unmap的比对
        # 第一种情况下，我们的bed文件是最全的，在确认part的时候我们的bothRep文件所构建的index还没有迁入新的序列，所以bed文件是包含所有的BEP index所有的序列
        # 第二种情况下，
        # 答案是 可以，首先bed文件是没有经过改动的那个文件，这里的left序列是在left 更新过的index里面，也就是UpIndex/UpIndex里面，这个index是根据
        # 几次的去掉序列而非增加序列所构建的leftRepclutot的那个文件来做的。所以left bed是包含所有的index的left文件
        RefChr = bedInfo[refContigName][0]
        RefStart = bedInfo[refContigName][1]
        RefEnd = bedInfo[refContigName][2]
        startJudge = int(RefStart) if int(RefStart)-2000>=0 else 0
        EndJudge = int(RefEnd) + 2000
        contigMapInfo = []
        # 这个也没有问题，我们的mateLinkFile是包含所有的序列，是最全的
        '''
        for lineMate in MateLinkFile.readlines():
            lineMate = lineMate.strip().split()
            if str(lineMate[1] + "-" + lineMate[len(lineMate)-1]) == queryContigName:
                contigMapInfo.append([lineMate[6],lineMate[8],lineMate[9]])
        '''
        log.info(f"range-for-ref{refContigName}-is-start-{RefStart}-end-{RefEnd}")

        for contig in MateLinkDirc:
            if contig == queryContigName:
                log.info(f"get-the-anchor-read-for-{queryContigName}-totally-" + str(len(MateLinkDirc[contig])) + "-alignments")
                for info in MateLinkDirc[contig]:
                    contigMapInfo.append(info)
        count = 0
        for i in contigMapInfo:
            if i[0] == RefChr and int(i[1]) >= startJudge and int(i[2]) <= EndJudge:
                count = count + 1
                log.info(f"suitble-anchor-found-{queryContigName}-chr-{i[0]}-start-{i[1]}-end-{i[2]}")
        if len(contigMapInfo) >= 5 and count >= (0.25 * len(contigMapInfo)) :
            outFile.write(str(refContigName + " " + queryContigName + "\n"))
    BedFile.close()
    outFile.close()
    MateLinkFile.close()
    candiFile.close()


