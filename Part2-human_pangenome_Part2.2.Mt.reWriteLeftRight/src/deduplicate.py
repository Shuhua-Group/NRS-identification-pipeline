import os
import collections
import networkx as nx
from Bio import SeqIO

def best_hit_alignment(alignment_path, identity, coverage, outdir):
    f_write = open(str(outdir + "/rmRedunt.tsv"), 'w')
    f = open(alignment_path, 'r')
    record = 2
    hit_contigs = set()
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')

        cov1 = 1.0 * (float(line[5]) - int(line[12])) / float(line[4])
        cov2 = 1.0 * (float(line[5]) - int(line[12])) / float(line[3])

        '''
        Blast file 
        1 2 3....
        a b
        a c
        a d
        a e 
        f q 
        f p
        f z
        .......等等
        '''

        contig1 = line[0].split("=")
        contig2 = line[1].split("=")

        if contig1[0] == contig1[1] and contig2[0] == contig2[1]:
            if line[0] not in hit_contigs:
                if line[0] == line[1]:
                    record = 1
                else:
                    record = 2
                    hit_contigs.add(line[0])
                    if float(line[2]) >= identity and (cov1 >= coverage and cov2 >= coverage):
                        for num in range(len(line) - 1):
                            f_write.write(line[num] + '\t')
                        f_write.write(line[len(line) - 1] + '\n')
            elif line[0] in hit_contigs and record == 1:
                record = 2
                hit_contigs.add(line[0])
                if float(line[2]) >= identity and (cov1 >= coverage and cov2 >= coverage):
                    for num in range(len(line) - 1):
                        f_write.write(line[num] + '\t')
                    f_write.write(line[len(line) - 1] + '\n')
    f_write.close()
    f.close()




def placement_pos(BEP_bed, LEP_bed, REP_bed, ArgOutDir):
    posUp2 = {}
    for record in SeqIO.parse(str(ArgOutDir + "/leftEnding/UpRepCluTot2.fa"), "fasta"):
        id = str(record.id).split("=")[0]
        posUp2[id] = "-l"
    for record in SeqIO.parse(str(ArgOutDir + "/rightEnding/UpRepCluTot2.fa"), "fasta"):
        id = str(record.id).split("=")[0]
        posUp2[id] = "-r"
    for record in SeqIO.parse(str(ArgOutDir + "/bothEnding/UpRepCluTot2.fa"), "fasta"):
        id = str(record.id).split("=")[0]
        posUp2[id] = "-b"

    f = open(BEP_bed, 'r')
    contig_pos = {}
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        name = line[3].split("-")[0] + "-" + line[3].split("-")[1]
        start = int(line[1])
        end = int(line[2])
        if name in posUp2:
            contig_pos[name] = [str(line[0]), start, end, posUp2[name]]
    # pos.add(name)
    f.close()

    f = open(LEP_bed, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        name = line[3].split("-")[0] + "-" + line[3].split("-")[1]
        start = int(line[1])
        end = int(line[2])
        if name in posUp2:
            contig_pos[name] = [str(line[0]), start, end, posUp2[name]]
    f.close()

    f = open(REP_bed, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        name = line[3].split("-")[0] + "-" + line[3].split("-")[1]
        start = int(line[1])
        end = int(line[2])
        if name in posUp2:
            contig_pos[name] = [str(line[0]), start, end, posUp2[name]]
    f.close()

    return contig_pos


def judge_distance(distance, contig_pos, pass_alignment, outdir):
    obtain_contigs = collections.defaultdict(list)
    map_contigs = []
    f = open(str(outdir + "/rmRedunt.tsv"), 'r')
    f_write = open(str(pass_alignment + ".filter"), 'w')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        q1 = line[0].split("=")[0]
        q2 = line[1].split("=")[0]

        if q1 not in contig_pos:
            name1 = q1.split('-')[0]
            name2 = q1.split('-')[1]
            query1_chr = contig_pos[name1][0]
            query1_begin = contig_pos[name1][1]
            query1_end = [name2][1]
            orentation1 = '-b'

        else:
            query1_chr = contig_pos[q1][0]
            query1_begin = contig_pos[q1][1]
            query1_end = contig_pos[q1][2]
            orentation1 = contig_pos[q1][3]

        if q2 not in contig_pos:
            name1 = q2.split('-')[0]
            name2 = q2.split('-')[1]
            query2_chr = contig_pos[name1][0]
            query2_begin = contig_pos[name1][1]
            query2_end = contig_pos[name2][1]
            orentation2 = '-b'
        else:
            query2_chr = contig_pos[q2][0]
            query2_begin = contig_pos[q2][1]
            query2_end = contig_pos[q2][2]
            orentation2 = contig_pos[q2][3]

        if query2_chr == query1_chr:
            if query2_begin - query1_end < distance and query1_begin - query2_end < distance:
                if q2 not in obtain_contigs:
                    map_contigs.append([q1 + orentation1, q2 + orentation2])
                    obtain_contigs[q1].append(q2)
                else:
                    if q1 not in obtain_contigs[q2]:
                        map_contigs.append([q1 + orentation2, q2 + orentation2])

                        obtain_contigs[q1].append(q2)

    f.close()
    G = nx.Graph()
    G.add_nodes_from(sum(map_contigs, []))
    info = [[(s[i], s[i + 1]) for i in range(len(s) - 1)] for s in map_contigs]
    for i in info:
        G.add_edges_from(i)
    for i in nx.connected_components(G):
        information = ''
        for z in i:
            information += z + '\t'
        f_write.write(information + '\n')  # for i in nx.connected_components(G)]
    f_write.close()


def rmRedundance(alignment_path, distance, identity, coverage, ArgOut):
    best_hit_alignment(alignment_path, identity, coverage, str(ArgOut + "/AllPlaced"))
    BEP_bed = ArgOut + "/bothEnding/Both.ending.bed"
    LEP_bed = ArgOut + "/leftEnding/Left.ending.bed"
    REP_bed = ArgOut + "/rightEnding/Right.ending.bed"
    contig_pos = placement_pos(BEP_bed, LEP_bed, REP_bed, ArgOut)
    judge_distance(distance, contig_pos, alignment_path, str(ArgOut + "/AllPlaced"))
