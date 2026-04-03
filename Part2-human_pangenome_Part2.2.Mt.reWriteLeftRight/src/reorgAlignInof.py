import os
import shutil

'''
'''
'''
def judge(file):
    left = {}
    right = {}
    f = open(file, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        identity = float(line[6])
        cov1 = float(line[9])
        cov2 = float(line[10])
        left_contig = line[11]
        right_contig = line[12]
        left_length = float(line[7])
        right_length = float(line[8])


        if left_contig not in left and right_contig not in left and left_contig not in right:
            if right_contig not in right:
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            else:
                if right[right_contig][2] < identity:
                    del left[right[right_contig][0]]

                    # left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,right_length]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                elif right[right_contig][2] == identity and right[right_contig][3] < cov1 and right[right_contig][
                    4] < cov2:
                    # left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,right_length]
                    del left[right[right_contig][0]]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
        elif left_contig in left and right_contig not in left and left_contig not in right:

            # right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            if left[left_contig][2] < identity:
                del right[left[left_contig][1]]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left[left_contig][2] == identity and left[left_contig][3] < cov1 and left[left_contig][4] < cov2:
                del right[left[left_contig][1]]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
        elif right_contig in left and left_contig not in right :
            if right_contig not in right:
                if right[left_contig][2] < identity or (left[left_contig][2] == identity and left[left_contig][3] < cov1 and left[left_contig][4] < cov2):
                    del right[left[right_contig][1]]
                    del left[right_contig]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
        elif left_contig in right and right_contig not in left:
            del left[right[left_contig][0]]
            del right[left_contig]
            left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
    f.close()
    return left, right
'''
def judge(file):
    left = {}
    right = {}
    f = open(file, 'r')
    contigPool = []
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        identity = float(line[6])
        cov1 = float(line[9])
        cov2 = float(line[10])
        left_contig = line[11]
        right_contig = line[12]
        left_length = float(line[7])
        right_length = float(line[8])

        if left_contig not in contigPool and right_contig not in contigPool:
            right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            contigPool.append(left_contig)
            contigPool.append(right_contig)
        elif left_contig in contigPool and right_contig not in contigPool:
            contigPool.append(right_contig)
            if left_contig in left:
                if left[left_contig][2] < identity or (left[left_contig][2] == identity and left[left_contig][3] < cov1 and left[left_contig][4] < cov2):
                    del right[left[left_contig][1]]

                    # left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,right_length]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left_contig in right:
                if right[left_contig][2] < identity or (right[left_contig][2] == identity and right[left_contig][3] < cov1 and right[left_contig][4] < cov2):
                    del left[right[left_contig][0]]
                    del right[left_contig]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
        elif right_contig in contigPool and left_contig not in contigPool:

            if right_contig in left:
                if left[right_contig][2] < identity or (left[right_contig][2] == identity and left[right_contig][3] < cov1 and left[right_contig][4] < cov2):
                    del right[left[right_contig][1]]
                    del left[right_contig]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif right_contig in right:
                if right[right_contig][2] < identity or (right[right_contig][2] == identity and right[right_contig][3] < cov1 and right[right_contig][4] < cov2):
                    del left[right[right_contig][0]]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
        elif right_contig in contigPool and left_contig in contigPool :
            if left_contig in left and right_contig in left:
                if (left[left_contig][2] < identity or (left[left_contig][2] == identity and left[left_contig][3] < cov1 and left[left_contig][4] < cov2)) and (left[right_contig][2] < identity or (left[right_contig][2] == identity and left[right_contig][3] < cov1 and left[right_contig][4] < cov2)) :
                    del right[left[left_contig][1]]
                    del right[left[right_contig][1]]
                    del left[left_contig]
                    del left[right_contig]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left_contig in right and right_contig in right:
                if (right[left_contig][2] < identity or (right[left_contig][2] == identity and right[left_contig][3] < cov1 and right[left_contig][4] < cov2)) and (right[right_contig][2] < identity or (right[right_contig][2] == identity and right[right_contig][3] < cov1 and right[right_contig][4] < cov2)) :
                    del left[right[left_contig][0]]
                    del left[right[right_contig][0]]
                    del right[left_contig]
                    del right[right_contig]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left_contig in left and right_contig in right:
                if (left[left_contig][2] < identity or (left[left_contig][2] == identity and left[left_contig][3] < cov1 and left[left_contig][4] < cov2)) and (right[right_contig][2] < identity or (right[right_contig][2] == identity and right[right_contig][3] < cov1 and right[right_contig][4] < cov2)):

                    del right[left[left_contig][1]]
                    del left[right[right_contig][0]]
                    del left[left_contig]
                    del right[right_contig]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left_contig in right and right_contig in left:
                if (right[left_contig][2] < identity or (right[left_contig][2] == identity and right[left_contig][3] < cov1 and right[left_contig][4] < cov2)) and (left[right_contig][2] < identity or (left[right_contig][2] == identity and left[right_contig][3] < cov1 and left[right_contig][4] < cov2)):
                    del left[right[left_contig][0]]
                    del right[left[right_contig][1]]
                    del left[right_contig]
                    del right[left_contig]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]

    f.close()
    return left, right


def alignment_priority(LEP_bed, REP_bed, path):
    #Identity_path, Contained_path, Overlap_path, Part_align_path = paths[0],paths[1],paths[2],paths[3]
    Identity_path, Contained_path, Overlap_path, Part_align_path = [str(path + "/identity"), str(path + "/contain"),str(path + "/overlap"), str(path + "/part")]
    upIde, upCon, upOver, upPart = [str(path + "/UPidentity"), str(path + "/UPcontain"),str(path + "/UPoverlap"), str(path + "/UPpart")]
    f = open(LEP_bed, 'r')
    reverse_info = {}
    pos_begin = {}
    pos_end = {}
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        seq_name = line[3].split("-")[0] + "-" +  line[3].split("-")[1]
        reverse = line[4]
        reverse_info[seq_name] = reverse
        pos_begin[seq_name] = line[1]
        pos_end[seq_name] = line[2]
    f.close()

    f = open(REP_bed, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        seq_name = line[3].split("-")[0] + "-" +  line[3].split("-")[1]
        reverse = line[4]
        reverse_info[seq_name] = reverse
        pos_begin[seq_name] = line[1]
        pos_end[seq_name] = line[2]
    f.close()
    identityPair = []

    f = open(Identity_path, 'r')
    left_identity = set()
    right_identity = set()
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        left_contig = line[11] 
        right_contig = line[12]
        left_identity.add(left_contig)
        right_identity.add(right_contig)
        identityPair.append((left_contig , right_contig))
    f.close()

    left_contain = set()
    right_contain = set()
    f_contain = open(upCon, 'w')
    f = open(Contained_path, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        left_contig = line[11]
        right_contig = line[12]
        if left_contig not in left_identity and right_contig not in right_identity and (left_contig, right_contig) not in identityPair:
        #if (left_contig, right_contig) not in identityPair:
            for num in range(len(line) - 1):
                f_contain.write(line[num] + '\t')
            f_contain.write(line[len(line) - 1] + '\n')
        left_contain.add(left_contig)
        right_contain.add(right_contig)
        identityPair.append((left_contig,right_contig))
    f_contain.close()
    f.close()
    left_identity = left_identity | left_contain
    right_identity = right_identity | right_contain


    left_overlap = set()
    right_overlap = set()
    f_overlap = open(upOver, 'w')
    f = open(Overlap_path, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        left_contig = line[11]
        right_contig = line[12]
        #if reverse_info[left_contig] == reverse_info[right_contig] and pos_begin[left_contig] <= pos_end[right_contig]: 
        if left_contig not in left_identity and right_contig not in right_identity and (left_contig, right_contig) not in identityPair:
        #if (left_contig, right_contig) not in identityPair:
            for num in range(len(line) - 1):
                f_overlap.write(line[num] + '\t')
            f_overlap.write(line[len(line) - 1] + '\n')
        left_overlap.add(left_contig)
        right_overlap.add(right_contig)
        identityPair.append((left_contig, right_contig))
    f_overlap.close()
    f.close()
    left_identity = left_identity | left_overlap
    right_identity = right_identity | right_overlap


    left_part = set()
    right_part = set()
    f_part = open(upPart, 'w')
    f = open(Part_align_path, 'r')
    for line in f.readlines():
        line = line.split("\n")[0]
        line = line.split("\t")
        left_contig = line[11]
        right_contig = line[12]
        if left_contig not in left_identity and right_contig not in right_identity and (left_contig, right_contig) not in identityPair:
        #if (left_contig, right_contig) not in identityPair:
            for num in range(len(line) - 1):
                f_part.write(line[num] + '\t')
            f_part.write(line[len(line) - 1] + '\n')
    f.close()
    f_part.close()

def redefine_alignment_result(Identity_path, save_folder):

    f_identity = open(save_folder+"final_identity.txt", 'w')
    left, right = judge(str(Identity_path))
    for key in left:
        f_identity.write(str(left[key][0]) + '\t' + str(left[key][1]) + '\t')
        if left[key][5] > left[key][6]:
            f_identity.write(str(left[key][0]) + '-l' + '\n')
        else:
            f_identity.write(str(left[key][1]) + '-r' + '\n')
    f_identity.close()



    f_contain = open(save_folder+"final_contained.txt", 'w')
    #f = open("/dev/shm/contained1.txt", 'r')
    left, right = judge(str(save_folder + "/UPcontain"))
    '''
        if left_contig not in left:
            if right_contig not in right:
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            else:
                
                #if right[right_contig][2] < identity:
                #    left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,
                #                                    right_length]
            
                
                if right[right_contig][2] < identity:
                    del left[right[right_contig][0]]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                elif right[right_contig][2] == identity and (right[right_contig][3] < cov1 and right[right_contig][4] < cov2):
                   #left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,right_length]
                    del left[right[right_contig][0]]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]

        else:
            #right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            if left[left_contig][2] < identity:
                del right[left[left_contig][1]]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]

                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left[left_contig][2] == identity and (left[left_contig][3] < cov1 and left[left_contig][4] and cov2):
                del right[left[left_contig][1]]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]

                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
    f.close()
    '''
    for key in left:
        f_contain.write(str(left[key][0]) + '\t' + str(left[key][1]) + '\t')
        if left[key][5] > left[key][6]:
            f_contain.write(str(left[key][0]) + '-l' + '\n')
        else:
            f_contain.write(str(left[key][1]) + '-r' + '\n')
    f_contain.close()

    f_overlap = open(save_folder +"final_overlap.txt", 'w')
    #f = open(str(save_folder + "/UPoverlap"), 'r')
    left, right = judge(str(save_folder + "/UPoverlap"))
    for key in left:
        f_overlap.write(str(left[key][0]) + '\t' + str(left[key][1]) + '\t')
        if left[key][5] > left[key][6]:
            f_overlap.write(str(left[key][0]) + '-l' + '\n')
        else:
            f_overlap.write(str(left[key][1]) + '-r' + '\n')
    f_overlap.close()




    f_part = open(save_folder +"final_part.txt", 'w')

    left, right = judge(str(save_folder + "/UPpart"))

    for key in left:
        f_part.write(str(left[key][0]) + '\t' + str(left[key][1]) + '\t')
        if left[key][5] > left[key][6]:
            f_part.write(str(left[key][0]) + '-l' + '\n')
        else:
            f_part.write(str(left[key][1]) + '-r' + '\n')
    f_part.close()

    '''
    if os.path.exists("/dev/shm/overlap1.txt"):
        os.remove("/dev/shm/overlap1.txt")
    if os.path.exists("/dev/shm/contained1.txt"):
        os.remove("/dev/shm/contained1.txt")
    if os.path.exists("/dev/shm/part1.txt"):
        os.remove("/dev/shm/part1.txt")
    '''
def conFirmPart(mergePath, config):
    IniEnv = str(config.getValue("IniEnv"))
    upPartFile = open(str(mergePath + "/confirm.final_part.txt"),'w')
    with open(str(mergePath + "/final_part.txt"),'r') as f:
        for line in f.readlines():
            confirmed = False
            line = line.strip().split()
            file = str(line[0] + "-" + line[1])
            os.system(str(IniEnv + " nucmer -p " + str(mergePath + "/" + file) + "/Lrep_Rcluster " + str(mergePath + "/" + file) + "/" + line[0] + ".fa" + " " + str(mergePath + "/" + file) + "/" + line[1] + ".cluster.fa" ))
            os.system(str(IniEnv + " nucmer -p " + str(mergePath + "/" + file) + "/Rrep_Lcluster " + str(mergePath + "/" + file) + "/" + line[1] + ".fa" + " " + str(mergePath + "/" + file) + "/" + line[0] + ".cluster.fa"))
            os.system(str(IniEnv + "delta-filter -r -q -g " + str(mergePath + "/" + file) + "/Lrep_Rcluster.delta > " + str(mergePath + "/" + file) + "/filter.Lrep_Rcluster.delta"))
            os.system(str(IniEnv + "delta-filter -r -q -g " + str(mergePath + "/" + file) + "/Rrep_Lcluster.delta > " + str(mergePath + "/" + file) + "/filter.Rrep_Lcluster.delta"))
            os.system(str(IniEnv + " show-coords -H -T -l -c -o " + str(mergePath + "/" + file) + "/filter.Lrep_Rcluster.delta > " + str(mergePath + "/" + file) + "/filter.Lrep_Rcluster.coord"))
            os.system(str(IniEnv + " show-coords -H -T -l -c -o " + str(mergePath + "/" + file) + "/filter.Rrep_Lcluster.delta > " + str(mergePath + "/" + file) + "/filter.Rrep_Lcluster.coord"))
            coordFile1 = open(str(mergePath + "/" + file + "/filter.Lrep_Rcluster.coord"), 'r')
            for Coordline1 in coordFile1:
                Coordline1 = Coordline1.strip("\n").split()
                if float(Coordline1[6]) > 97:
                    confirmed = True
            coordFile2 = open(str(mergePath + "/" + file + "/filter.Rrep_Lcluster.coord"), 'r')
            for Coordline2 in coordFile2:
                Coordline2 = Coordline2.strip("\n").split()
                if float(Coordline2[6]) > 97:
                    confirmed = True
            if confirmed == True:
                upPartFile.write('\t'.join(line) + "\n")




def PopinsMerge(mergePath, config):
    if not os.path.isdir(mergePath + "/popinsMerge"):os.mkdir(mergePath + "/popinsMerge")
    finalOverlap = open(str(mergePath + "/final_overlap.txt"),'r')
    for line in finalOverlap:
        line = line.strip()
        # left - right
        line = line.split()
        outdir = str(mergePath + "/popinsMerge/" + line[0] + "-" + line[1])
        if not os.path.isdir(outdir):os.mkdir(outdir)
        if not os.path.isdir(outdir + "/" + line[0]): os.mkdir(outdir + "/" + line[0])
        if not os.path.isdir(outdir + "/" + line[1]): os.mkdir(outdir + "/" + line[1])
        shutil.copyfile(str(mergePath + "/" + line[0] + "-" + line[1] + "/" + line[0] + ".fa"), str(outdir + "/" + line[0] + "/contigs.fa"))
        shutil.copyfile(str(mergePath + "/" + line[0] + "-" + line[1] + "/" + line[1] + ".fa"), str(outdir + "/" + line[1] + "/contigs.fa"))
        Inienv = "module load bifrost sickle samtools/1.9 bwakit popins/2.0.0;"

        popinsMergeRes = os.system(str(Inienv + "popins merge -p " + outdir  + " -c " + outdir + "/" + line[0] + "-" + line[1] + ".fa"))
        if not popinsMergeRes == 0:
            print("popins Merge is not successful; break down, you need to get the UpRepCluTot2.fa yourself ")

    finalContain = open(str(mergePath + "/final_contained.txt"), 'r')
    for line in finalContain:
        line = line.strip()
        # left - right
        line = line.split("\t")
        outdir = str(mergePath + "/popinsMerge/" + line[0] + "-" + line[1])
        if not os.path.isdir(outdir):os.mkdir(outdir)
        if not os.path.isdir(outdir + "/" + line[0]): os.mkdir(outdir + "/" + line[0])
        if not os.path.isdir(outdir + "/" + line[1]): os.mkdir(outdir + "/" + line[1])
        shutil.copyfile(str(mergePath + "/" + line[0] + "-" + line[1] + "/" + line[0] + ".fa"), str(outdir + "/" + line[0] + "/contigs.fa"))
        shutil.copyfile(str(mergePath + "/" + line[0] + "-" + line[1] + "/" + line[1] + ".fa"), str(outdir + "/" + line[1] + "/contigs.fa"))
        Inienv = "module load bifrost sickle samtools/1.9 bwakit popins/2.0.0;"
        #Inienv = "source /anaconda3/bin/activate ; conda activate panGenomePart2 ;"
        popinsMergeRes = os.system(str(Inienv + "popins merge -p " + outdir  + " -c " + outdir + "/" + line[0] + "-" + line[1] + ".fa"))
        if not popinsMergeRes == 0:
            print("popins Merge is not successful; break down, you need to get the UpRepCluTot2.fa yourself ")
