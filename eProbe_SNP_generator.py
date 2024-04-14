import argparse
import logging
from eModules.pipe_controler import *
from eModules.fasta_operator import *
import os

usage = """python eProbe_SNP_generator.py -d SNPs_df.txt -g genome.fasta -l 81 -s shift -o output"""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument("-f", "--snp_df", dest="snp_df", action="store", nargs='?',
                    help="Input SNPs data frame (tsv).",
                    metavar="FILE", required=True)
parser.add_argument("-r", "--reference", dest="reference", action="store", nargs='?',
                    help="The reference genome (fasta) used to generate probes. "
                         "The same one used for VCF calling.",
                    metavar="FILE", required=True)
parser.add_argument("-t", "--thread", dest="thread", action="store", nargs='?',
                    help="Number of thread (default: 1).",
                    default="1", metavar="STRING")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output prefix for process files and result files.",
                    metavar="STRING", required=True)
parser.add_argument("--method", dest="method", action="store", nargs='?',
                    help="Method for extracting seq from genome (pos/edge, default: pos)."
                         "pos: base on the position of SNPs."
                         "edge: base on the start and end of testing probes (ignore all other parameter).",
                    choices=['pos', 'edge'], default="pos", metavar="STRING")
parser.add_argument("--probe_length", dest="probe_length", action="store", nargs='?',
                    help="Length of probes (odd number recommended).",
                    default="81", metavar="STRING")
parser.add_argument("--probe_shift", dest="probe_shift", action="store", nargs='?',
                    help="Position shift (default: 0) fo method pos. - for shifting left "
                         "and + for shifting right.",
                    default="0", metavar="STRING")
parser.add_argument("--replace", dest="replace", action="store", nargs='?',
                    help="Replace target SNP with a base neither ref not alt (on/off, default: off).",
                    choices=['on', 'off'], default="off", metavar="STRING")
parser.add_argument("--type", dest="type", action="store", nargs='?',
                    help="Select probes based on mutation type (ts/tv, default: both).",
                    choices=["ts", "tv"], default=None, metavar="STRING")
args = parser.parse_args()

SNP_DF = args.snp_df
REFERENCE = args.reference
THREAD = args.thread
OUTPUT_PREFIX = args.output
METHOD = args.method
PROBE_LENGTH = args.probe_length
PROBE_SHIFT = args.probe_shift
REPLACE = args.replace
TYPE = args.type

# make script searching list
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "eModules")
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

# make log file
logging.basicConfig(filename=f"{OUTPUT_PREFIX}.generator.log", level=logging.INFO)

generator_script_name = "Get_probes_from_SNPs.py"
generator_params = ["-d", SNP_DF, "-g", REFERENCE, "-l", PROBE_LENGTH, "-s", PROBE_SHIFT, "-m", METHOD, "--replace",
                    REPLACE, "--type", TYPE, "-o", OUTPUT_PREFIX]
run_python_script(generator_script_name, PATH_SEARCH_LIST, generator_params)
