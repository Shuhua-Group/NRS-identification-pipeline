import os
import collections
import networkx as nx
from Bio import SeqIO
# 对这个脚本做一些改变，原来的含义不变, 或者就按照原来的意思来做也没什么问题，只是一些地方需要改动一下。
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
        # 下面这一段代码，就是，我们如果成会一直读取，直到找到某一个序列的可以匹配的第二个序列，
        # 比如，如果我们的blast文件是
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
        # 我们会读取 a b 然后判断a b是否是我们想要的，如果是我们想要的，那么对于a开头的就不再往下读了，就直接把a 写到我们的hit_contigs
        # 如果a b不符合条件就过，然后继续往下找。
        # 然后把我们的这个blast的结果写到f_write里面，算是做了一个筛选/
        contig1 = line[0].split("=")
        contig2 = line[1].split("=")
        # 代表rep
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


# 下面这里我们把position全部都读进去也反映了一点，那就是我们的align的操作是针对所有的contig的，left；right以及both的。

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
        # 这一种情况对应的是 merge的情况，也就是我们的序列是被merge的
        # 这里假想的情况是我们如果left和right需要merge的话，我们会记我们的representative为 left(xx=xx)-right(xx==xx)
        # 但是我在实际的操作的时候是直接 是选择了left或者right中最长的一个作为我们的代表序列。比较简单，而且实际上的效果应该是笔他这里用的好
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
# 整个流程的运行来看，基本的想法就是我们把序列相互比对，然后比的上的，而且在2000bp内的represent序列就归纳到一个集合里面来。


def rmRedundance(alignment_path, distance, identity, coverage, ArgOut):
    best_hit_alignment(alignment_path, identity, coverage, str(ArgOut + "/AllPlaced"))
    BEP_bed = ArgOut + "/bothEnding/Both.ending.bed"
    LEP_bed = ArgOut + "/leftEnding/Left.ending.bed"
    REP_bed = ArgOut + "/rightEnding/Right.ending.bed"
    contig_pos = placement_pos(BEP_bed, LEP_bed, REP_bed, ArgOut)
    judge_distance(distance, contig_pos, alignment_path, str(ArgOut + "/AllPlaced"))
