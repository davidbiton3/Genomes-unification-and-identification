import os
import sys
import linecache
from datetime import datetime
from package import alg


# this function will read a file by given path + filename.
# # return the file data.
def read_file(file_name):
    
    data = []
    try:
        with open(file_name) as file:
            for line in file:
                data.append(line.split())
    except FileNotFoundError:
        print_error_to_log(True, "")
        print("cant open file, exit from the program. ")
        sys.exit(0)

    return data


# this function run all the alg with his needed args.
# didnt return a value, just run.
def run(data_input, result_folder_path, n_threads, completeness_threshold, contamination_threshold,
        genome_group_size, ani_threshold, linkage_threshold):
    # now = datetime.now()
    # date_time = now.strftime("%d-%m-%Y_%H%M%S")
    # final_path = result_folder_path + "/result_" + str(date_time)
    # try:
    #     if not os.path.isdir(final_path):
    #         os.mkdir(final_path)
    #     os.chdir(final_path)  # change directory
    # except OSError:
    #     print_error_to_log(True, "")
    #     print("please check your path: " + result_folder_path)
    #     sys.exit(1)

    final_path = result_folder_path + "/result_23-11-2020_211126"
    if os.path.isdir(final_path):
        os.chdir(final_path)  # change directory
        print("its change dir from func run")
    alg.alg(data_input, n_threads, completeness_threshold, contamination_threshold,
            genome_group_size, ani_threshold, linkage_threshold)


# this function write to the log file all the error.
def print_error_to_log(is_exception, message):
    now = datetime.now()
    date_time = "[" + str(now.strftime("%d-%m-%Y %H:%M")) + "] "

    log_file = open("log.txt", 'a')
    if is_exception:
        exc_type, exc_obj, tb = sys.exc_info()
        f = tb.tb_frame
        line_no = tb.tb_lineno
        filename = f.f_code.co_filename
        filename = filename.split("/")[-1]
        linecache.checkcache(filename)
        info_string = str(date_time) + (
            ' ERROR FROM EXCEPTION: (file: {}, line: {}): {}'.format(filename, line_no, exc_obj))
        log_file.write(info_string + '\n')
    else:
        log_file.write(str(date_time) + message + '\n')

    log_file.close()


# this function write to the log file all the info.
def print_info_to_log(message):
    now = datetime.now()
    date_time = "[" + str(now.strftime("%d-%m-%Y %H:%M")) + "] INFO:\t"
    log_file = open("log.txt", 'a')
    log_file.write(str(date_time) + message + '\n')
    log_file.close()


# this function will check if path is a file and file isn't empty.
# return a boolean as answer.
def is_non_zero_file(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


# this function will check line from the input file of the program
# return boolean and a message.
def check_is_legal_input_row(input_row):
    
    n = len(input_row)
    if n != 3:
        return False, "legal row is 3 columns only!"

    if input_row[1][-4] == ".faa" or input_row[2][-4] == ".fna":
        return False, "legal row is tree columns like (genome_name tab path.fna tab path.faa)"

    return True, ""
