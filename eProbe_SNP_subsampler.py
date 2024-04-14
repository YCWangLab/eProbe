import argparse
from eModules.pipe_controler import *
from eModules.fasta_operator import *
import os
import logging
usage = """python eProbe_SNP_subsampler.py -f snp.tsv -r ref.fasta -t threads -o output"""
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-f", "--snp_df", dest="snp_df", action="store", nargs='?',
                    help="Input SNPs data frame (tsv).",
                    metavar="FILE", required=True)
parser.add_argument("-r", "--reference", dest="reference", action="store", nargs='?',
                    help="The reference genome (fasta) used to generate probes. "
                         "The same one used for VCF calling.",
                    metavar="FILE", required=True)
parser.add_argument("--window_size", dest="window_size", action="store", nargs='?',
                    help="Window size for sub-sampling (default: 1000). "
                         "eProbe will select one SNP from each non-overlapping window containing SNPs.",
                    default="1000", metavar="STRING")
parser.add_argument("-t", "--thread", dest="thread", action="store", nargs='?',
                    help="Number of thread (default: 1).",
                    default="1", metavar="STRING")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output prefix for process files and result files.",
                    metavar="STRING", required=True)

parser.add_argument("--desire_Tm", dest="desire_tm", action="store", nargs='?',
                    help="Desire melting temperature (default: 75).",
                    default="75", metavar="STRING")
parser.add_argument("--desire_GC", dest="desire_gc", action="store", nargs='?',
                    help="Desire GC content (default: 50).",
                    default="50", metavar="STRING")
parser.add_argument("--select_weights", dest="weights", action="store", nargs='?',
                    help="Weights for Tm, GC, hairpin, dimer, DUST (default: 0.5,0.5,0,0,0)."
                         "eProbe will first calculate the score according to the distance of each"
                         "biophysical property to their desired value (for hairpin, dimer and DUST score, "
                         "the smaller the better), then select the SNP with the best overall performance."
                         "If not provided, eProbe will randomly select one SNP from a window.",
                    metavar="STRING")

parser.add_argument("--keep_bed", dest="keep_bed", action="store", nargs='?',
                    help="Annotation file. Usually a 4-column BED file converted from gff, "
                         "with last column describing the function of this genomic region. "
                         "If provided, a SNP of a window will be kept prioritarily if it overlaps with BED.",
                    metavar="FILE")
parser.add_argument("--keep_bed_column", dest="keep_bed_column", action="store", nargs='?',
                    help="Number of the BED file column used to match feature. "
                         "If not provided, keep all genomic regions of the BED file.",
                    metavar="STRING")
parser.add_argument("--keep_bed_feature", dest="keep_bed_feature", action="store", nargs='?',
                    help="What feature in the column want to be kept (e.g., exon). "
                         "If not provided, keep all genomic regions of the BED file.",
                    metavar="STRING")

parser.add_argument("--output_mode", dest="output_mode", action="store", nargs='?',
                    help="Output mode (both/feature/window, default: window)."
                         "window: If SNPs in a window falling within the BED (if provided), prioritize selecting one "
                         "from these SNPs. Otherwise, then either select one based on weights (if provided), or "
                         "randomly choose one."
                         "feature: Only output all SNPs falling within the BED."
                         "both: Output both Window and Feature SNP sets and their merged SNP set",
                    choices=["both", "feature", "window"], default="window", metavar="STRING")

parser.add_argument("--probe_number", dest="probe_number", action="store", nargs='?',
                    help="Randomly subsample SNPs to a certain number."
                         "If the selected SNPs number is less than the input SNPs number, output selected SNPs."
                         "If output mode is Both, subsample based on merged SNPs dataframe.",
                    default=None, metavar="STRING")
parser.add_argument("--seed", dest="seed", action="store", nargs='?',
                    help="Seed number for random subsampling.",
                    default="123", metavar="STRING")

args = parser.parse_args()
# common input:
SNP_DF = args.snp_df
REFERENCE = args.reference
OUTPUT_PREFIX = args.output
THREAD = args.thread

DESIRE_TM = args.desire_tm
DESIRE_GC = args.desire_gc
SELECT_WEIGHTS = args.weights
WINDOW_SIZE = args.window_size
KEEP_BED = args.keep_bed
KEEP_BED_COLUMN = args.keep_bed_column
KEEP_BED_FEATURE = args.keep_bed_feature
SELECT_MODE = args.output_mode

SEED = args.seed
PROBE_NUMBER = args.probe_number
SUBSAMPLE_DF1 = args.subsample_df1
SUBSAMPLE_DF2 = args.subsample_df2

# make log file
logging.basicConfig(filename=f"{OUTPUT_PREFIX}.subsampler.log", level=logging.INFO)

# make script searching list
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "eModules")
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

# generate chr_len.tsv
generate_chr_length_tsv(read_fasta(REFERENCE), OUTPUT_PREFIX)

# run SNPs_selector.py
selector_script_name = "SNPs_selector.py"
selector_params = ["-d", SNP_DF, "-t", DESIRE_TM, "-g", DESIRE_GC, "-w", SELECT_WEIGHTS, "-c",
                   f"{OUTPUT_PREFIX}_chr_len.tsv", "-s", WINDOW_SIZE, "-b", KEEP_BED, "-n", KEEP_BED_COLUMN, "-f",
                   KEEP_BED_FEATURE, "-m", SELECT_MODE, "-r", THREAD, "-o", OUTPUT_PREFIX]
run_python_script(selector_script_name, PATH_SEARCH_LIST, selector_params)

if PROBE_NUMBER:
    # run SNPs_subsampler.py
    subsampler_script_name = "SNPs_subsampler.py"
    if SELECT_MODE == "window":
        subsampler_params = ["-s", SEED, "-r", PROBE_NUMBER, "--df1", f"{OUTPUT_PREFIX}_window_SNPs.tsv", "-o",
                             OUTPUT_PREFIX]
    elif SELECT_MODE == "feature":
        subsampler_params = ["-s", SEED, "-r", PROBE_NUMBER, "--df1", f"{OUTPUT_PREFIX}_feature_SNPs.tsv", "-o",
                             OUTPUT_PREFIX]
    elif SELECT_MODE == "both":
        subsampler_params = ["-s", SEED, "-r", PROBE_NUMBER, "--df1", f"{OUTPUT_PREFIX}_merged_SNPs.tsv", "-o",
                             OUTPUT_PREFIX]
    else:
        raise ValueError("Input proper select mode: both, feature, window")
