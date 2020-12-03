import os
import argparse
import general_functions


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    # construct the argument parse and parse the arguments
    ap.add_argument("-i", "--inputFile", required=True,
                    help=" input file - mandatory flag ")
    ap.add_argument("-p", "--pathRoute", required=True,
                    help=" result folder path -  mandatory flag ")
    ap.add_argument("-n", "--nThreads", type=int, default=20, required=False,
                    help=" number of cpus - optional flag ")
    ap.add_argument("-comp", "--completeness", type=float, default=85.0, required=False,
                    help=" completeness threshold - optional flag")
    ap.add_argument("-cont", "--contamination", type=float, default=5.0, required=False,
                    help=" contamination threshold - optional flag ")
    ap.add_argument("-s", "--genomeSize", type=int, default=5000, required=False,
                    help=" genome group size - optional flag ")
    ap.add_argument("-a", "--ani", type=float, default=95.0, required=False,
                    help=" ani threshold - optional flag ")
    ap.add_argument("-l", "--linkage", type=float, default=70.0, required=False,
                    help=" linkage threshold - optional flag ")
    args = vars(ap.parse_args())

    input_file_name = args["inputFile"]
    result_folder_path = args["pathRoute"]
    n_threads = int(args["nThreads"])
    completeness_threshold = float(args["completeness"])
    contamination_threshold = float(args["contamination"])
    genome_group_size = int(args["genomeSize"])
    ani_threshold = float(args["ani"])
    linkage_threshold = float(args["linkage"])

    if os.path.exists("log.txt"):
        os.remove("log.txt")

    data_input = general_functions.read_file(input_file_name)
    general_functions.run(data_input, result_folder_path, n_threads, completeness_threshold,
                          contamination_threshold, genome_group_size, ani_threshold, linkage_threshold)
