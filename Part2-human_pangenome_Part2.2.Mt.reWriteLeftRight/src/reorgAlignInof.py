import os
import shutil

'''
这个函数做了什么呢？
首先这个函数读取了contig的一些信息，这个信息是前两步我们的分析的结果就是我们输出的那个文件，就是我写的那个流程的分析结果，这个结果包含了contig的
在染色体上insert的信息。
然后我们依次的读取left 和right map的结果，然后把map的结果一边写到我们的left identity； right identity两个set中去，一边把我们的结果给输出出来到文件中去
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

        # 下面是我复制粘贴的上面的代码， 再下面''' '''内的是源代码
        if left_contig not in left and right_contig not in left and left_contig not in right:
            if right_contig not in right:
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            else:
                if right[right_contig][2] < identity:
                    del left[right[right_contig][0]]
                    # 把已经有的left的东西替换掉。这里给了替换的两个条件，第一个是identity，第二个是identity结合cov
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
            # 同上；也是判断以及替换，最终是选择最合适的一对pair，比如 left1;left2都能和right1匹配，那么可能要判断left1或者2里面那个能比较好的和right1匹配，然后选择他作为代表
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
        # 下面是我复制粘贴的上面的代码， 再下面''' '''内的是源代码
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
                    # 把已经有的left的东西替换掉。这里给了替换的两个条件，第一个是identity，第二个是identity结合cov
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
            # 同上；也是判断以及替换，最终是选择最合适的一对pair，比如 left1;left2都能和right1匹配，那么可能要判断left1或者2里面那个能比较好的和right1匹配，然后选择他作为代表
            # right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
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
                    # 这里不会出现同一个left和right在同一个identity中出现的情况，因为之前的设定，第一次写identity等文件的时候是不允许这种情况的
                    # 所以这里不会出现把right中相关的删除以后，再用right来找找不到的情况。
                    # 但是事实就是之前没有过滤好
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
        left_contig = line[11] # 这里仍然拿到的是left和right的比对的结果，show-coord的结果
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

    # 这里还是会有这个问题，我们能不能通过bed文件来判断我们的比对的rep的位置，是可以的，rep还是那句话，只少不多，所以可以直接用。
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
                
                #上面是源代码，我改动了一下，如果现在的pair的情况比较好，我们就舍弃掉原来的pair
                #例如：如果原来是 l1-r1 , 我们的新的pair是 l2 - r1，如果l2和r1的匹配更好，那么我们就直接把 r1匹配的那一个写成l2，然后把
                # left里面 l1和r1的key-value pair删掉，然后同时在left的字典里面写上 l2-r1
                
                if right[right_contig][2] < identity:
                    del left[right[right_contig][0]]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                elif right[right_contig][2] == identity and (right[right_contig][3] < cov1 and right[right_contig][4] < cov2):
                   #left[right[right_contig][0]] = [left_contig, right_contig, identity, cov1, cov2, left_length,right_length]
                    del left[right[right_contig][0]]
                    right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                    left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
        # 下面的代码这样写的原因是我们在同一个类比如identity里面没有l1-r1-match1; l1-r1-match2这种情况，所以，如果left相同，right必定不同
        else:
            #right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            if left[left_contig][2] < identity:
                del right[left[left_contig][1]]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                #right的改动
                right[right_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
            elif left[left_contig][2] == identity and (left[left_contig][3] < cov1 and left[left_contig][4] and cov2):
                del right[left[left_contig][1]]
                left[left_contig] = [left_contig, right_contig, identity, cov1, cov2, left_length, right_length]
                #right的改动
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
        #Inienv = "source /Users/kongshuang/opt/anaconda3/bin/activate ; conda activate panGenomePart2 ;"
        popinsMergeRes = os.system(str(Inienv + "popins merge -p " + outdir  + " -c " + outdir + "/" + line[0] + "-" + line[1] + ".fa"))
        if not popinsMergeRes == 0:
            print("popins Merge is not successful; break down, you need to get the UpRepCluTot2.fa yourself ")
