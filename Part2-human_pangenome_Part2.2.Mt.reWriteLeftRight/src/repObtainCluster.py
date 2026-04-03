import argparse
import os

from Bio import SeqIO
import logAndcon



def rep_cluster(seq_path, path_merge_bed, path_rep, path_contig, fullContigPath):
    seq = {}

    for record in SeqIO.parse(seq_path, 'fasta'):
        seq[str(record.id)] = [str(record.seq), str(record.seq.reverse_complement())]
    f = open(path_merge_bed, 'r')
    fullContigFile = open(fullContigPath, 'w')
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        length = int(line[2]) - int(line[1])
        contigs = line[3]
        contigs = contigs.split(',') 
        contig_name1 = []
        for contig_info in contigs:
            name = contig_info.split("-")
            contig_name = name[0]
            sample_name = name[1]
            '''
            for num in range(2, len(name) - 1):
                sample_name += name[num]
                sample_name += '_'
            sample_name += name[len(name) - 1]  
            '''                                    
            #reverse1 = contig_info.split(".")[1].split('_')[1]
            #length1 = contig_info.split(".")[1].split('_')[0]
            reverse1 = contig_info.split("-")[2]
            length1 = contig_info.split("-")[3]
            if str(contig_info.split("-")[0] + "-" + contig_info.split("-")[1]) in seq:
                contig_name1.append([int(length1), str(contig_info.split("-")[0] + "-" + contig_info.split("-")[1]), contig_name, sample_name, reverse1])
        contig_name = sorted(contig_name1, key=(lambda x: x[0]), reverse=True) 
        if len(contig_name) != 0:
            R_name = contig_name[0][1]
            reverse = contig_name[0][4]
            f_write_R = open(path_rep + R_name + ".fa", 'w')

            if reverse == "0":
                f_write_R.write('>' + R_name + '\n')
                f_write_R.write(seq[R_name][0] + '\n')
            elif reverse == "1":
                f_write_R.write('>' + R_name + '\n')
                f_write_R.write(seq[R_name][1] + '\n')
            f_write_R.close()

            f_write_contig = open(path_contig + R_name + '.fa', 'w')
            for contig_num in range(len(contig_name)):
                if contig_name[contig_num][4] == "0":
                    f_write_contig.write('>' + contig_name[contig_num][1] + '\n')
                    f_write_contig.write(seq[contig_name[contig_num][1]][0] + '\n')

                    fullContigFile.write('>' + contig_name[contig_num][1] + '\n')
                    fullContigFile.write(seq[contig_name[contig_num][1]][0] + '\n')

                elif   contig_name[contig_num][4] == "1":
                    f_write_contig.write('>' + contig_name[contig_num][1]+ '\n')
                    f_write_contig.write(seq[contig_name[contig_num][1]][1] + '\n')

                    fullContigFile.write('>' + contig_name[contig_num][1]+ '\n')
                    fullContigFile.write(seq[contig_name[contig_num][1]][1] + '\n')
            f_write_contig.close()
    f.close()

def nucmerRep(repFile, clusterPath, IniEnv):
    nucPair = {}
    repFile = repFile[:-1]
    slash = repFile.rfind('/',2)
    rootPath = repFile[0:slash]
    RepCompareDir = rootPath + "/RepCompare"
    if not os.path.isdir(RepCompareDir): os.mkdir(RepCompareDir)
    upCluster = rootPath + "/UpRepCluster"
    if not os.path.isdir(upCluster): os.mkdir(upCluster)
    repTotal = rootPath + "/RepTotal"
    UpRepClusterTotal = rootPath + "/UpRepCluTot"
    UpRepClusterTotalFile = open(UpRepClusterTotal, 'w')

    repTotalFile = open(repTotal, 'w')
    for root, dirs, files in os.walk(repFile):
        for file in files:
            # test: nucPair[str(root + file)] = str(clusterFile + file)
            #print(nucPair[str(root + file)])
            repSeq = open(str(root + "/" + file), 'r')
            for line in repSeq.readlines():
                repTotalFile.write(line)
            RepFile = file[0:(len(file)-3)]
            ''' 
            **************************************************
            **************************************************
            TODO : don't forget to redo this command line 
            **************************************************
            **************************************************
            '''
            os.system(str(IniEnv + "nucmer -p " + RepCompareDir + "/" + RepFile + " "  + str(root + "/" + file + " ") + str(clusterPath + file)))
            #**************************************************TESTTESTTESTTESTTESTTESTTEST
            #os.system(IniEnv + "show-coords -H -T -l -c -o " + RepCompareDir + "/" + RepFile + " > " +  RepCompareDir + "/" +  RepFile + ".coord")
            deltaFile = open(str(RepCompareDir + "/" + RepFile + ".delta"), 'r')
            clusterFile = open(str(clusterPath + file), 'r')
            clusterRead = []
            upClusterFile = upCluster + "/" + file
            upClusterFileFa = open(upClusterFile,'w')
            for line in deltaFile.readlines():
                if line[0] == '>':
                    line = line.strip()
                    clusterRead.append(line.split()[1])
            header = 0
            for line in clusterFile.readlines():
                if line[0] == ">":
                    line = line.strip()
                    if line[1:len(line)] in clusterRead:
                        header = 1
                        upClusterFileFa.write(line + "\n")
                        UpRepClusterTotalFile.write(line + "=" + RepFile + "\n")
                    else:
                        header = 0

                else:
                    if header == 1 :
                        upClusterFileFa.write(line)
                        UpRepClusterTotalFile.write(line)
                        header = 0 