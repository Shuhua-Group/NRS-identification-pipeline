import re
import logging
log = logging.getLogger(__name__)

#src.clusterCompare.EnsureCandidate(str(bothEndDir + "/bothLeftCandidate.txt"), Args.outdir, str(bothEndDir + "/Both.ending.bed"), str(bothEndDir + "/UpBothLeftCandi"))
def EnsureCandidate(CandiFile, rootPath, bedFile, outfile): 
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
 
        RefChr = bedInfo[refContigName][0]
        RefStart = bedInfo[refContigName][1]
        RefEnd = bedInfo[refContigName][2]
        startJudge = int(RefStart) if int(RefStart)-2000>=0 else 0
        EndJudge = int(RefEnd) + 2000
        contigMapInfo = []

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


