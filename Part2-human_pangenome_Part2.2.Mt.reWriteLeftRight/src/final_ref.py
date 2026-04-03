import os
from Bio import SeqIO
from collections import Counter
import collections


def overlap_contigs(placed_align_path, ArgOutdir):
    #contig_len = {}
    contig_len = {}
    with open(str(ArgOutdir + "/leftEnding/Left.ending.bed"), 'r') as f:
        for line in f.readlines():
            line = line.split('\n')[0]
            line = line.split('\t')
            contigInfo = line[3].split("-")
            contig_len[str(contigInfo[0] + "-" + contigInfo[1])] = int(contigInfo[3])
    f.close()

    with open(str(ArgOutdir + "/rightEnding/Right.ending.bed"), 'r') as f:
        for line in f.readlines():
            line = line.split('\n')[0]
            line = line.split('\t')
            contigInfo = line[3].split("-")
            contig_len[str(contigInfo[0] + "-" + contigInfo[1])] = int(contigInfo[3])
    f.close()

    with open(str(ArgOutdir + "/bothEnding/Both.ending.bed"), 'r') as f:
        for line in f.readlines():
            line = line.split('\n')[0]
            line = line.split('\t')
            contigInfo = line[3].split("-")
            contig_len[str(contigInfo[0] + "-" + contigInfo[1])] = int(contigInfo[3])
    f.close()
    contigCluster = []
    contig_cluster = {}
    contig_type = {}
    delete_contigs = set()
    f = open(placed_align_path, 'r')

    for line in f.readlines():
        cluster = set()
        line = line.split('\n')[0]
        line = line.strip()
        line = line.split('\t')
        length = 0
        label = 'all'

        for contig in line:

            # if len(contig) > 2:
            if len(line) >=2 :
                ID = contig.split('-')[0] + "-" + contig.split('-')[1]
                if contig_len[ID] >= length:
                    if label != 'b':
                        length = contig_len[ID]
                        rep = ID
                        type = contig.split('-')[-1]
                        if type == 'b':
                            label = 'b'
                    else:
                        if contig.split('-')[-1] == 'b':
                            length = contig_len[ID]
                            rep = ID
                cluster.add(ID)
                delete_contigs.add(ID)
        contig_cluster[rep] = cluster
        contig_type[rep] = type
    f.close()
    #    print len(contig_cluster)
    return contig_cluster, contig_type, delete_contigs,


def obtain_contigs(rep ,cluster ,RepCluTot):
    # cluster and rep
    for record in SeqIO.parse(RepCluTot, "fasta"):
        id = str(record.id).split("=")
        if id[0] == id[1]:
            rep[id[1] + "=" + id[1]].append(str(record.seq))
            cluster[id[1] + "=" + id[1]].append([str(record.id),str(record.seq)])
        else:
            cluster[id[1] + "=" + id[1]].append([str(record.id),str(record.seq)])
    return rep, cluster


def generate_placed_contigs(contig_cluster, contig_type, delete_contigs, left_rep, right_rep, both_rep, cluster,leftFinalUpRepCluTot, rightFinalUpRepCluTot, bothFinalUpRepCluTot):


    right_rep_path = open(rightFinalUpRepCluTot, 'w')
    left_rep_path = open(leftFinalUpRepCluTot, 'w')
    both_rep_path = open(bothFinalUpRepCluTot, 'w')
    
	
    for key in right_rep:
        if key.split("=")[1] not in delete_contigs:
            for subcontig in (cluster[key]):
                right_rep_path.write('>' + subcontig[0]  + '\n')
                right_rep_path.write(subcontig[1] + '\n')


    for key in left_rep:
        if key.split("=")[1]  not in delete_contigs:
            for subcontig in cluster[key]:
                left_rep_path.write('>' + subcontig[0] + '\n')
                left_rep_path.write(subcontig[1] + '\n')


    for key in both_rep:
        if key.split("=")[1]  not in delete_contigs:
            for subcontig in cluster[key]:
                both_rep_path.write('>' + subcontig[0] + '\n')
                both_rep_path.write(subcontig[1] + '\n')


    for contig in contig_cluster:
        #        print contig_type[contig]
        if contig_type[contig] == 'b':
            seq = both_rep[str(contig + "=" + contig)][0]
   #         both_rep_path.write('>' + contig + '\n')
    #        both_rep_path.write(seq + '\n')

            for contig_name in (contig_cluster[contig]):
                for subcontig in cluster[str(contig_name + "=" + contig_name)]:
                    both_rep_path.write('>' + subcontig[0] + '\n')
                    both_rep_path.write(subcontig[1] + '\n')
        elif contig_type[contig] == 'l':
            seq = left_rep[str(contig + "=" + contig)][0]
    #        left_rep_path.write('>' + contig + '\n')
     #       left_rep_path.write(seq + '\n')
            for contig_name in (contig_cluster[contig]):
                for subcontig in cluster[str(contig_name + "=" + contig_name)]:
                    left_rep_path.write('>' + subcontig[0] + '\n')
                    left_rep_path.write(subcontig[1] + '\n')
        elif contig_type[contig] == 'r':
            seq = right_rep[contig][0]
      #      right_rep_path.write('>' + contig + '\n')
       #     right_rep_path.write(seq + '\n')
            for contig_name in (contig_cluster[contig]):
                for subcontig in cluster[str(contig_name + "=" + contig_name)]:
                    right_rep_path.write('>' + subcontig[0] + '\n')
                    right_rep_path.write(subcontig[1] + '\n')

    right_rep_path.close()
    left_rep_path.close()
    right_rep_path.close()


def finalRef(ArgOutDir):
    left_rep = collections.defaultdict(list)

    cluster = collections.defaultdict(list)
    right_rep = collections.defaultdict(list)
    both_rep = collections.defaultdict(list)
    placed_align_path = str(ArgOutDir + "/AllPlaced/AllPlacedAligned.tsv.filter")
    contig_cluster, contig_type, delete_contigs = overlap_contigs(placed_align_path, ArgOutDir)
    leftRepCluTot = ArgOutDir + "/leftEnding/UpRepCluTot2.fa"
    rightRepCluTot = ArgOutDir + "/rightEnding/UpRepCluTot2.fa"
    bothRepCluTot = ArgOutDir + "/bothEnding/UpRepCluTot2.fa"

    left_rep,cluster = obtain_contigs(left_rep, cluster, leftRepCluTot)
    right_rep,cluster = obtain_contigs(right_rep, cluster, rightRepCluTot)
    both_rep,cluster = obtain_contigs(both_rep, cluster, bothRepCluTot)
    leftFinalUpRepCluTot = str(ArgOutDir + "/leftEnding/FinalUpRepCluTot")
    rightFinalUpRepCluTot = str(ArgOutDir + "/rightEnding/FinalUpRepCluTot")
    bothFinalUpRepCluTot = str(ArgOutDir + "/bothEnding/FinalUpRepCluTot")
    generate_placed_contigs(contig_cluster, contig_type, delete_contigs, left_rep, right_rep, both_rep, cluster,leftFinalUpRepCluTot, rightFinalUpRepCluTot,bothFinalUpRepCluTot)


#def rmRedundancyUnpla():
    
