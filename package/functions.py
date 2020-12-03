import os
import sys
import itertools
import general_functions as function
from package.Genome import Genome
from package.Cluster import Cluster


# this function will check the quality of the genomes.
# didnt return a value, but will create an out file called "checkm.stdout".
def check_m(data_input, n_threads):

    function.print_info_to_log("# state 1 #")
    try:
        if os.path.isdir('checkm'):
            os.rmdir("checkm")

        os.mkdir("checkm")
        os.mkdir("checkm/proteins")

        # create soft-link to all proteins files
        for index, protein in enumerate(data_input):
            is_valid, error_message = function.check_is_legal_input_row(protein)
            if not is_valid:
                message = "ERROR FROM IF:(file: functions.py, line: 23): input row isn't valid " \
                          "(row index in file: "+str(index) + ") " + error_message
                function.print_error_to_log(False, message)
                print("input file error, check log file.")
                sys.exit(0)

            path = protein[2]
            name_gen = protein[0] + ".faa"
            os.system("ln -s " + path + " checkm/proteins/" + name_gen)

    except OSError:
        function.print_error_to_log(True, "")

    function.print_info_to_log("starting checkm")
    os.system("~itaish/software/bin/checkm.lineage_wf.sh checkm/proteins checkm --ncpus=" + str(n_threads))
    function.print_info_to_log("ending checkm")


# this function will initialize the genomes form the given data arg.
# didnt return a value.
def init_genomes(data_input, completeness_threshold, contamination_threshold):
    function.print_info_to_log("# state 2 #")
    if os.path.isfile('checkm/checkm.stdout'):
        with open('checkm/checkm.stdout.txt', 'w') as outfile:
            with open("checkm/checkm.stdout", 'r') as infile:
                outfile.write(infile.read())
    else:
        message = "ERROR FROM IF:(file: functions.py, line: 47): file not exist "
        function.print_error_to_log(False, message)
        print("can't open file, exit from the program.")
        sys.exit(0)

    file_genome = function.read_file("checkm/checkm.stdout.txt")[3:-1]
    start_len_file = len(file_genome)
    if not start_len_file:
        message = "ERROR FROM IF:(file: functions.py, line: 55): file is empty "
        function.print_error_to_log(False, message)
        print("file is empty, exit from the program.")
        sys.exit(0)
    function.print_info_to_log("starting the mapping genomes process")
    genome_list = []
    for index, line in enumerate(file_genome):
        completeness = float(line[12])
        contamination = float(line[13])
        if completeness >= completeness_threshold and contamination <= contamination_threshold:
            name = data_input[index][0]
            path = data_input[index][1]
            genome_list.append(Genome(name, path, completeness, contamination))

    removed_genomes = start_len_file - len(genome_list)
    function.print_info_to_log("ending the mapping genomes process - " + str(removed_genomes) + " genomes was removed")

    return genome_list


# this function will sorted the genomes list by the genome quality
# in decreasing order.
# didnt return a value, just sort.
def sorted_genome(genome_list):
    function.print_info_to_log("# state 3 #")
    function.print_info_to_log("# sorted the genomes by quality, decreasing")
    sorted(genome_list, key=lambda t: t.quality, reverse=True)
    return genome_list


# this function will clustering all the given genomes_list to clusters.
# will cluster by taking a subgroup (from top) and run fastANI against the the bigger group the include the subgroup.
# return list of all created clusters.
def cluster_the_genomes(genome_list, genome_group_size, ani_threshold, n_threads):
    function.print_info_to_log("# state 4 #")
    # helper list of tuples, let us know if genome is clustered
    genome_list_bool = list(map(lambda gen: (gen.path, False), genome_list))
    clusters_list = []
    genome_to_remove = []
    need_to_cluster = True
    # state a
    if not os.path.isdir("init_clusters"):
        os.mkdir("init_clusters")

    index = 0
    function.print_info_to_log("starting to clustering the genomes")
    while need_to_cluster:

        path = "init_clusters/round-" + str(index)
        index += 1  # up the index to the around iteration
        if not os.path.isdir(path):
            os.mkdir(path)  # create dir by iteration.

        if len(genome_list) < genome_group_size:
            genome_group_size = len(genome_list)

        queries_list = list(itertools.islice(genome_list, genome_group_size))

        create_path_file_from_genome(path + "/queries.list.txt", queries_list, False)

        # state b
        create_path_file_from_genome(path + "/subjects.list.txt", genome_list, False)

        # state c
        outfile = path + "/fastANIOutFile.txt"
        # os.system(
        #     "fastANI --ql " + (path + "/queries.list.txt") + " --rl " + (path + "/subjects.list.txt") + " -o " +
        #     outfile + " -t " + str(n_threads) + " > " + outfile + ".stderr 2> " + outfile + ".stdout")
        outfile = function.read_file(outfile)

        if not outfile[-1]:  # remove last empty line from the output file.
            outfile = outfile[:-1]

        if outfile:  # file not empty.
            # state d
            running_outfile = True
            while running_outfile:  # this loop is running until the file is clear.
                left_genome = outfile[0][0]

                if left_genome == outfile[0][1]:
                    outfile = outfile[1:]

                # partition  all left_genome only  and right_genome.
                sub_file, outfile = partition(lambda sub_row: (sub_row[0] == left_genome) and not (
                        sub_row[0] == left_genome and sub_row[1] == left_genome), outfile)

                if not outfile:
                    running_outfile = False

                # if is_cluster_left return false its not clustered.
                is_cluster_left, index_left = check_if_genome_clustered(left_genome, genome_list_bool)

                if not is_cluster_left:
                    genome = genome_list[index_left]  # get the genome.
                    genome_list_bool[index_left] = (genome_list_bool[index_left][0], True)
                    if index_left >= genome_group_size:
                        genome_to_remove.append(index_left)
                    temp_cluster = Cluster(genome)  # create single cluster.

                    if sub_file:  # sub_file - is all left_genome exiting.
                        for left_row_file in sub_file:
                            right_genome = left_row_file[1]
                            is_cluster_right, index_right = check_if_genome_clustered(right_genome, genome_list_bool)
                            if not is_cluster_right:  # if right path is cluster if not can be cluster.
                                ani_score = float(left_row_file[2])
                                if ani_score >= ani_threshold:
                                    temp_cluster.genome_list.append(genome_list[index_right])
                                    genome_list_bool[index_right] = (genome_list_bool[index_right][0], True)
                                    if index_right >= genome_group_size:
                                        genome_to_remove.append(index_right)
                    else:
                        function.print_info_to_log("sub_file is empty, the genome name: "+genome.name+" is a singleton")
                    clusters_list.append(temp_cluster)
        else:
            message = "ERROR FROM IF:(file: functions.py, line: 131): file is empty"
            function.print_error_to_log(False, message)
            print("file is empty, exit from the program. ")
            sys.exit(0)

        # state e
        genome_list = update_the_list(genome_list, genome_to_remove, genome_group_size)
        genome_list_bool = update_the_list(genome_list_bool, genome_to_remove, genome_group_size)
        genome_to_remove = []
        need_to_cluster = len(genome_list) > 0

    function.print_info_to_log("ending to clustering the genomes")
    return clusters_list


# this function calculate the link score of all genomes and mapping those genomes that didnt compliance to conditions.
# will create a new folder to each cluster, then will run fastANI all X all,
# func: calc_link_of_cluster - will calculate the link score of every genome of clusters,
# func: maping_clusters - will mapping all those genome that didnt compliance to conditions.
def calc_and_maping(clusters_list, n_threads, ani_threshold, linkage_threshold):
    function.print_info_to_log("# state 5 #")
    try:
        if not os.path.isdir("clusters"):
            os.mkdir("clusters")
    except OSError:
        function.print_error_to_log(True, "")

    cluster_singleton_list_to_add = []
    function.print_info_to_log("# starting calculate link and maping")
    for cluster in clusters_list:
        len_genome_list = len(cluster.genome_list)
        path = "clusters/cluster_" + str(cluster.id)
        try:
            if not os.path.isdir(path):
                os.mkdir(path)
        except OSError:
            function.print_error_to_log(True, "")

        genome_list = cluster.genome_list

        all_path_left_genome = []
        name_file = path + "/genome_path_list_before_maping.txt"
        create_path_file_from_genome(name_file, genome_list, False)

        if len_genome_list > 1:

            outfile = path + "/fastANIOutFile.txt"
            # run fastANi Against himself
            # os.system("fastANI --ql " + name_file + " --rl " + name_file + " -o " + outfile + " -t "
            #           + str(n_threads) + " > " + outfile + ".stderr 2> " + outfile + ".stdout")

            if not init_out_file(outfile, len_genome_list, all_path_left_genome):
                continue

            #  starting calc
            calc_score_linkage_list_genome(cluster.genome_list, all_path_left_genome, ani_threshold, path)
            for x in cluster.genome_list:
                if x.linked <= 0:
                    print("true"+str(cluster.id))
                    print(len(cluster.genome_list))

            # create a friendship file before the mapping state.
            write_friendship_outfile(path + "/friendship_before_maping.txt", cluster.genome_list)

            # maping the linked score of all the genomes in the cluster less the linkage_threshold (is_low_linked).
            temp_singleton_cluster_list_to_add = maping_cluster(cluster, linkage_threshold)
            for x in cluster.genome_list:
                if x.linked <= 0:
                    print("trueeeeeeeeee" + str(cluster.id))
                    print(len(cluster.genome_list))

            if len(cluster.genome_list) == 1:
                print("xxx" + str(genome_list[0].linked))
                print(cluster.id)
                print(len(temp_singleton_cluster_list_to_add))

            # if have change in on the cluster, add new singleton cluster to the list ,and create file to print.
            if temp_singleton_cluster_list_to_add:
                cluster_singleton_list_to_add += temp_singleton_cluster_list_to_add
                if len(cluster.genome_list) == 1:
                    print(len(temp_singleton_cluster_list_to_add))
                create_folder_to_singleton_cluster(temp_singleton_cluster_list_to_add, cluster.id)
                write_friendship_outfile(path + "/friendship_after_maping.txt", cluster.genome_list)

        create_path_file_from_genome(path + "/genome_path_list_after_maping.txt", cluster.genome_list, False)
        if cluster.id == 1026:
            return

    function.print_info_to_log("# ending calculate link clusters and maping - " +
                               str(len(cluster_singleton_list_to_add)) + " new singleton cluster was created")

    if cluster_singleton_list_to_add:
        clusters_list += cluster_singleton_list_to_add

    return clusters_list


# this function will select a representative to the cluster,
# by called method of class "Cluster", select_representative()
# return a list of all representatives that selected.
def select_reperesentive(clusters_list):
    function.print_info_to_log("# state 6 #")
    function.print_info_to_log("starting selecting clusters representatives")
    representatives_list = []
    for cluster in clusters_list:
        if len(cluster.genome_list) > 1:
            cluster.select_reperesentive()

        representatives_list.append(cluster.representative_genome)

    function.print_info_to_log("ending selecting clusters representatives")

    return representatives_list


# this function will preform a gtdbtk state, that helper to classification the
# belonging of the representative genome to the tree of life
# didnt return a value, just set the cluster the taxonomy of the representative genome.
def gtdbtk(clusters_list, representatives_list):
    function.print_info_to_log("# state 7 #")
    function.print_info_to_log("create representatives_list file")

    if not os.path.isdir("gtdbtk"):
        os.mkdir("gtdbtk")

    create_path_file_from_genome("gtdbtk/representatives_list.txt", representatives_list, True)

    out_dir = "gtdbtk"
    function.print_info_to_log("running the gtdbtk comment")

    os.system("gtdbtk classify_wf --batchfile gtdbtk/representatives_list.txt --out_dir gtdbtk --cpus 20" + " > " +
              "gtdbtk/gtdbtk.stderr 2> " + "gtdbtk/gtdbtk.stdout")

    if (not os.path.isfile(out_dir + "/gtdbtk.ar122.summary.tsv")) and\
            (not os.path.isfile(out_dir + "/gtdbtk.bac120.summary.tsv")):
        message = "ERROR FROM IF:(file: functions.py): gtdbtk output files not exist!"
        function.print_error_to_log(False, message)
        print("gtdbtk output files not exist!, exit from the program. ")
        sys.exit(0)

    if os.path.isfile(out_dir + "/gtdbtk.ar122.summary.tsv"):
        outfile_one = function.read_file(out_dir + "/gtdbtk.ar122.summary.tsv")[1:]
        function.print_info_to_log("starting parse the the gtdbtk.ar122.summary.tsv file and create out file one.")
        for row_file_one in outfile_one:
            index = int(row_file_one[0])
            clusters_list[index].taxonomy = row_file_one[1]
        function.print_info_to_log("ending parse the the gtdbtk.ar122.summary.tsv file")

    if os.path.isfile(out_dir + "/gtdbtk.bac120.summary.tsv"):
        outfile_two = function.read_file(out_dir + "/gtdbtk.bac120.summary.tsv")[1:]
        function.print_info_to_log("starting parse the the gtdbtk.bac120.summary.tsv file")
        for row_file_two in outfile_two:
            index = int(row_file_two[0])
            clusters_list[index].taxonomy = row_file_two[1]
        function.print_info_to_log("ending parse the the gtdbtk.bac120.summary.tsv file")

    with open("representatives_taxonomy_list", 'w') as output:
        for cluster in clusters_list:
            output.write(str(cluster.id) + "\t" + cluster.representative_genome.name + "\t" + cluster.taxonomy + '\n')


# this function will create the first script output file.
# create a file with <genome-name>	<cluster-ID> for all genomes at all clusters.
# didnt return a value, just create a file.
def create_final_file(clusters_list):
    function.print_info_to_log("# state 8 #")
    with open("genome_cluster_list", 'w') as output:
        for cluster in clusters_list:
            for genome in cluster.genome_list:
                output.write(genome.name + "\t" + str(cluster.id) + "\n")

    function.print_info_to_log("finish...")


# this function will update the (list) by get a list of list index to remove.
# update_the_list - helper list of tuples, let us know if genome is clustered
# return a updated the list  without the tuples that all ready used.
def update_the_list(the_list, list_index_to_remove, top_group_size):
    the_list = the_list[top_group_size:]
    for index in list_index_to_remove:
        index = index - top_group_size
        del the_list[index]
    return the_list


# this function check genome is cluster
# return tuple (True/False ,index genome).
def check_if_genome_clustered(genome_path, genome_list_bool):
    try:
        index = genome_list_bool.index((genome_path, False))
        return False, index
    except ValueError:
        return True, -1


# this function will do partition at the given iterable, by the condition at the predicate.
# those element that compliance to the conditions stored at the trues list, else at the falses list.
# return trues and falses lists.
def partition(pred, iterable):
    trues = []
    falses = []
    for item in iterable:
        if pred(item):
            trues.append(item)
        else:
            falses.append(item)
    return trues, falses


# this function will create a friendship outfile file that explanation all the friendships inside the cluster genomes.
# didnt return a value, just create a file.friendship_outfile.
def write_friendship_outfile(file_name, genome_list):
    genome_list_len = len(genome_list)
    with open(file_name, 'w') as output:
        output.write("path \t linked \t percent link \n")
        for index, genome in enumerate(genome_list):
            genome_linked = genome.linked
            percent_link = str(0)
            if genome_list_len - 1:
                percent_link = str((genome_linked / (genome_list_len - 1)) * 100)
            output.write(genome.path + '\t' + str(genome_linked) + '\t' + percent_link)
            if index < genome_list_len - 1:
                output.write("\n")


# this function will get genome path and search the path in genome list
# return the index if found else -1.
def get_indexpath_from_genomes(item_path, genome_list):
    for index, gen in enumerate(genome_list):
        if item_path == gen.path:
            return index
    return -1


# this function will check if opposite is exist
# example have (a,b) check if exists (b,a) if yes return ani score of (b,a)
# all_path_left_genome - list of tuples, every tuple  will have on the left side (genome path) and
# right side (sub file that all is left side == genome path)
# return tuple (True/False, ani_score).
def bool_exists_opposite(left_path, right_path, all_path_left_genome):
    for tuple_element in all_path_left_genome:
        path_left = tuple_element[0]
        if path_left == right_path:
            for row in tuple_element[1]:
                if row[1] == left_path:
                    return True, row[2]
            return False, -1

    return False, -1


# this function will decrease friends from genome if not complain the condition
# didn't return a value, just update list friends on the genome.
def decrease_linkage(index_path_genome_linked, genome_path, genome_list):
    for index in index_path_genome_linked:
        if genome_list[index].delete_path_from_list_path(genome_path):
            genome_list[index].linked -= 1


# this function calc linkage and sava in list of tuple all friends with our ani score
# return the calc linkage and list of tuple all friends.
def calc_score_linkage_and_get_path_similar(ani_threshold, outfile):
    sum_link = 0
    list_genome_path = []
    for row in outfile:
        if float(row[2]) >= ani_threshold:
            sum_link += 1
            list_genome_path.append((row[1], float(row[2])))
    return sum_link, list_genome_path


# this function will mapping the genomes inside the cluster,
# that need to be delete because a marked as weak link.
# return cluster list after mapping.
def maping_cluster(cluster, linkage_threshold):
    cluster_list_to_add = []
    index_set_genome_kill = set()
    len_genome_list = len(cluster.genome_list)
    genome_list = cluster.genome_list

    need_to_map = True
    while need_to_map:
        need_to_map = False
        for j, genome in enumerate(genome_list):
            if j in index_set_genome_kill:
                continue
            if len_genome_list > 1:
                if genome.linked > 0:
                    if (genome.linked * 100) / (len_genome_list - 1) < linkage_threshold:
                        list_path_genome_linked, temp = zip(*genome.list_path_genome_linked)
                        index_path_genome_linked = get_indexs_from_list_path(list_path_genome_linked, genome_list)
                        decrease_linkage(index_path_genome_linked, genome.path, genome_list)
                        genome.update_linkage_list_path_genome_linked_is_low_linked(0, [])
                        index_set_genome_kill.add(j)
                        len_genome_list = len(genome_list) - len(index_set_genome_kill)
                        if len_genome_list > 1:
                            if len(index_path_genome_linked) > 0:
                                for index_value in index_path_genome_linked:  # check if friends is low linked
                                    linked = genome_list[index_value].linked
                                    if (linked * 100) / (len_genome_list - 1) < linkage_threshold:
                                        if j > index_value:
                                            need_to_map = True
                                            break

                        elif len_genome_list == 1:
                            index_remaining = int((set(range(0, len(cluster.genome_list))) -
                                                   index_set_genome_kill).pop())
                            cluster.genome_list[index_remaining].update_linkage_list_path_genome_linked_is_low_linked(
                                0, [])
                            index_set_genome_kill.add(index_remaining)
                            break

                        else:
                            break

                else:  # the genome isn't popular at the cluster, will remove.
                    genome.update_linkage_list_path_genome_linked_is_low_linked(0, [])
                    index_set_genome_kill.add(j)
                    len_genome_list = len(genome_list) - len(index_set_genome_kill)

            else:  # singleton cluster.
                index_set_genome_kill.add(j)
                genome.update_linkage_list_path_genome_linked_is_low_linked(0, [])
                break

    if len(index_set_genome_kill):  # genome will remove from the cluster.
        index_set_genome_kill = sorted(index_set_genome_kill, reverse=True)
        if len(index_set_genome_kill) == len(genome_list):
            del index_set_genome_kill[-1]  # this state save the original genome that create this cluster.

        for index in index_set_genome_kill:  # delete & create a singleton clusters to all genome killed.
            cluster_list_to_add.append(Cluster(genome_list[index]))
            del cluster.genome_list[index]
    # for x in cluster_list_to_add:
    #     print(x.id)
    return cluster_list_to_add


# this function will get lists path genome and return lists index genome
def get_indexs_from_list_path(list_path_genome_linked, genome_list):
    index_list_genome = []
    for index, path in enumerate(list_path_genome_linked):
        for genome in genome_list:
            if genome.path == path:
                index_list_genome.append(index)
                break

    return index_list_genome


# this function will create new file that contain only paths form genomes list.
# didnt return a value, create new file named "file_name".
def create_path_file_from_genome(file_name, genome_list, with_index):
    genome_list_len = len(genome_list)
    if not genome_list_len:
        message = "ERROR FROM IF:(file: functions.py, line: 502): genome list is empty"
        function.print_error_to_log(False, message)
        sys.exit(0)

    with open(file_name, 'w') as output:
        for index, genome in enumerate(genome_list):
            if with_index:
                output.write(genome.path + '\t' + str(index))
            else:
                output.write(genome.path)

            if index < genome_list_len - 1:
                output.write("\n")


# this function will initialize the "all path left genome" to all the genomes in the cluster.
# "all path left genome" is build from tuple (genome.path, all the couples with that genome from current outfile)
# return True if initialize was successful.
def init_out_file(outfile, len_genome_list, all_path_left_genome):
    outfile = function.read_file(outfile)
    if outfile:
        if not outfile[-1]:
            outfile = outfile[:-1]

        right_sub_file = []
        for index in range(len_genome_list):
            path_left = ""
            if index == 0:
                path_left = outfile[0][0]
                left_sub_file, right_sub_file = partition(
                    lambda sub_row: sub_row[0] == path_left and sub_row[1] != path_left, outfile)

            else:
                path_left = right_sub_file[0][0]
                left_sub_file, right_sub_file = partition(
                    lambda sub_row: sub_row[0] == path_left and sub_row[1] != path_left,
                    right_sub_file)

            # only if not exist state (A, A) at the top of right_sub_file.
            if right_sub_file[0][0] == right_sub_file[0][1]:
                right_sub_file = right_sub_file[1:]

            # create list of tuples, every tuple  will have on the left side (genome path) and
            # right side (sub file that all is left side == genome path)
            all_path_left_genome.append((path_left, left_sub_file))
    else:
        return False

    return True


# this function will calculate the score of genome in the cluster X all genomes at the cluster the if them will called
# friends, calc friendships by a given ani_score by the file outfile, (from fastANI)
# if the average score between genome A to genome B is bigger than the ani_threshold is count as a friend.
# print an "genome_not_have_opposite.txt" file that indicate a missing couple.
# didnt return a value, just calc.
def calc_score_linkage_list_genome(genome_list, all_path_left_genome, ani_threshold, path):

    with open(path + "/genome_not_have_opposite.txt", 'w') as outfile:
        for j, genome in enumerate(genome_list):
            genome_link, list_path_genome_linked = calc_score_linkage_and_get_path_similar(
                ani_threshold, all_path_left_genome[j][1])
            index_path_to_remove = []

            if len(list_path_genome_linked) > 0:
                for i, item_path in enumerate(list_path_genome_linked):  # check if have opposite.
                    is_opposite, ani_score_right = bool_exists_opposite(genome.path, item_path[0],
                                                                        all_path_left_genome)
                    if not is_opposite:
                        ani_score_right = item_path[1]
                        outfile.write(genome.path + '\t' + item_path[0] + '\n')

                    ani_score_left = item_path[1]
                    avg_score = (float(ani_score_left) + float(ani_score_right)) / 2
                    if avg_score < ani_threshold:  # need to be out me from friend.
                        genome_link -= 1  # out from me the friend.
                        index_path_to_remove.append(i)

                if index_path_to_remove:
                    index_path_to_remove.reverse()
                    for i_value in index_path_to_remove:
                        del list_path_genome_linked[i_value]

            genome.update_linkage_list_path_genome_linked_is_low_linked(genome_link, list_path_genome_linked)

    if not function.is_non_zero_file(path + "/genome_not_have_opposite.txt"):  # delete file if empty.
        os.remove(path + "/genome_not_have_opposite.txt")


# this function will create a folder for every singleton cluster.
# didnt return a value
def create_folder_to_singleton_cluster(singleton_cluster_list, cluster_id):
    print(str(cluster_id) + "   - "+str(len(singleton_cluster_list)))
    for cluster in singleton_cluster_list:
        path = "clusters/cluster_" + str(cluster.id)
        try:
            if not os.path.isdir(path):
                os.mkdir(path)
        except OSError:
            function.print_error_to_log(True, "")

        genome_list = cluster.genome_list

        name_file = path + "/genome_path.txt"
        create_path_file_from_genome(name_file, genome_list, False)
