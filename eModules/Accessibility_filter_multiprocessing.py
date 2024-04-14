# -*- coding:utf-8 -*-
from Bio import SeqIO
import os
import subprocess
from multiprocessing import Pool
import datetime
import shutil
from pipe_controler import *
from fasta_operator import *


def run_script(script_path, parameters):
    command = ['python', script_path] + parameters
    subprocess.check_output(command, universal_newlines=True)


def run_scripts_parallel(script_path, parameters_list, num_processes):
    with Pool(processes=num_processes) as pool:
        pool.starmap(run_script, [(script_path, params) for params in parameters_list])


def merge_fasta_files(fasta_paths, output_file):
    with open(output_file, 'w') as out_fasta:
        for fasta_path in fasta_paths:
            with open(fasta_path, 'r') as in_fasta:
                for line in in_fasta:
                    out_fasta.write(line)


if __name__ == '__main__':
    import argparse

    usage = """python Run_accessibility_filtering_multiprocessing.py -f probe.fasta -i bowtie2_index -t thread -o output"""

    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument("-f", "--fasta", dest="fasta", action="store", nargs='?',
                        help="Input fasta.", metavar="FILE", required=True)
    parser.add_argument("-i", "--index", dest="index", action="store", nargs='?',
                        help="Bowtie2 reference index.", metavar="FILE", required=True)
    parser.add_argument("-k", dest="keep_hit", action="store", nargs='?', default='100',
                        help="The number of hits to be kept in bowtie2 mapping (default: 100).", metavar="STRING")
    parser.add_argument("-t", "--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Thread for running bowtie2 mapping.", metavar="STRING")
    parser.add_argument("--mm", dest="m_mode", action="store", nargs='?', default="sensitive",
                        help="Bowtie2 local mapping mode (very-fast, fast, sensitive, very-sensitive). (default: sensitive)",
                        metavar="STRING")
    parser.add_argument("--trim5", dest="trim5", action="store", nargs='?', default="0",
                        help="Trim <int> bases from 5'/left end of reads (default: 0).",
                        metavar="STRING")
    parser.add_argument("--trim3", dest="trim3", action="store", nargs='?', default="0",
                        help="Trim <int> bases from 3'/right end of reads (default: 0).",
                        metavar="STRING")
    parser.add_argument("--fm", dest="f_mode", action="store", nargs='?', default="strict",
                        help="Accessibility filtering mode (strict, moderate, custom) (default: strict)."
                             "If custom is enabled, custom filtering thresholds will be used.",
                        metavar="STRING")
    parser.add_argument("--mapq", dest="mapq", action="store", nargs='?', default="30",
                        help="Minimum mapping quality (default: 30). Only effective in custom mode.", metavar="STRING")
    parser.add_argument("--nm", dest="nm", action="store", nargs='?', default="3",
                        help="Maximum number of edit distances (default: 3). Only effective in custom mode.",
                        metavar="STRING")
    parser.add_argument("--xo", dest="xo", action="store", nargs='?', default="0",
                        help="Maximum gap number (default: 0). Only effective in custom mode.", metavar="STRING")
    parser.add_argument("--xg", dest="xg", action="store", nargs='?', default="0",
                        help="Maximum total gap length (default: 0). Only effective in custom mode.", metavar="STRING")
    parser.add_argument("--xsr", dest="xsr", action="store", nargs='?', default="0.8",
                        help="Maximum secondary alignment score ratio (default: 0.8). Only effective in custom mode.",
                        metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="FILE", required=True)

    args = parser.parse_args()

    # making parameters_list
    # # command template: parameters that don't change in each subprocess
    temp_input = args.fasta
    temp_output = f"{args.output}.accessibility_filter_tempout"
    output_file = f"{args.output}.processing_probe.fasta"
    for db in args.index.split(","):
        # split the file and generate a temp_dir to store them
        print("Splitting probe file ...", datetime.datetime.now())
        split_files = split_fasta(temp_input, int(args.thread), args.output)
        print("Splitting probe file completed ...", datetime.datetime.now())
        # split_file = args.output + .split1.fasta

        parameter_template = []
        # # -t = 1 for each processing
        # # strict and moderate modes parameters: -i -k -t --mm --trim5 --trim3 --fm --nm --xo --xg --xsr
        if args.f_mode in ["strict", "moderate"]:
            parameter_template = ["-i", db, "-k", args.keep_hit, "--mm", args.m_mode, "--trim5", args.trim5,
                                  "--trim3", args.trim3, "--fm", args.f_mode]
        # # custom mode parameters: -k -t --mm --trim5 --trim3 --fm --mapq --nm --xo --xg --xsr
        else:
            parameter_template = ["-i", db, "-k", args.keep_hit, "--mm", args.m_mode, "--trim5", args.trim5,
                                  "--trim3", args.trim3, "--fm", args.f_mode, "--mapq", args.mapq, "--nm", args.nm,
                                  "--xo",
                                  args.xo, "--xg", args.xg, "--xsr", args.xsr]
        parameters_list = []
        # # multiprocessing parameters: -f -o
        # # make result_files
        result_files = []
        for split_file in split_files:
            prefix = split_file.replace(".fasta", "")
            result_file = f"{prefix}.temp_probe.fasta"
            result_files.append(result_file)
            parameters_list.append(
                parameter_template + ["-f", split_file, "-o", prefix])
        # 3 run python script multiprocessing
        script_directory = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.join(script_directory, 'Accessibility_probe_filter.py')
        print(f"Running accessibility filtering for db: {db} ...", datetime.datetime.now())
        run_scripts_parallel(script_path, parameters_list, int(args.thread))
        # 4 process results
        merge_fasta_files(result_files, temp_output)
        temp_input = f"{args.output}.accessibility_filter_tempin"
        os.rename(temp_output, temp_input)

    os.rename(temp_input, output_file)
    print("Accessibility filtering completed ...", datetime.datetime.now())
    before_filtering = str(int(int(subprocess.check_output(["wc", "-l", args.fasta]).split()[0]) / 2))
    after_filtering = str(int(diff_file_line_count(args.fasta, output_file) / 2))
    print("Filtered out %s probes from %s probes in accessibility filtering." % (after_filtering, before_filtering),
          datetime.datetime.now())

    # 5 remove temp dir
    shutil.rmtree(f"{args.output}_eProbe_temp")
