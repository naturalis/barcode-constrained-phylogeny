import os
import csv
import glob

def one_time_runners(files):
    one_time_dict = {}
    for filename in files:
        file = open(filename, "r").readlines()
        seconds = file[1].split("\t")[0]
        one_time_dict[filename.split("/")[1].replace(".benchmark.txt", "")] = seconds
    return one_time_dict


def parallel_runners(folders, time_dict):
    for folder in folders:
        benchmarks = glob.glob(os.path.join(folder, "*.benchmark.txt")) 
        for filename in benchmarks:
            file = open(filename, "r").readlines()
            seconds = file[1].split("\t")[0]
            time_dict[filename.split("/")[2].replace(".benchmark.txt", "")] = seconds
    return time_dict


def get_sequence_stats(filename):
    sequences_num = len(open(filename).readlines())/2
    n_count = 0
    total_count = 0
    with open(filename, 'r') as f:
        for count, line in enumerate(f, start=1):
            if count % 2 == 0:
                n_count = line.count("N")
                total_count = len(line)
    if not total_count == 0:
        n_percent = (n_count/total_count)*100
    else:
        n_percent = 0
    return sequences_num, n_percent

if __name__ == '__main__':
    path = "benchmarks/"
    benchmarks = glob.glob(os.path.join(path, "*.benchmark.txt")) 
    time_dict = one_time_runners(benchmarks)
    
    subfolders = [f.path for f in os.scandir(path) if f.is_dir()]
    time_dict = parallel_runners(subfolders, time_dict)
    time_dict = dict(sorted(time_dict.items()))

    with open('benchmarks/times.tsv', 'w') as csv_file:  
        writer = csv.writer(csv_file, delimiter='\t')
        for key, value in time_dict.items():
            if "-of-" in key:
                numbs = key.split("-")
                folder = f"{numbs[1]}-of-{numbs[3]}"
                sequences_num, n_percent = get_sequence_stats(f"results/fasta/family/{folder}/unaligned.fa")
                key = key.split("-")
                writer.writerow([key[0], key[1], value, int(sequences_num), n_percent])
            else:
                writer.writerow([key, value])
