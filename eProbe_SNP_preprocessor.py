import os
import subprocess
import argparse
from eModules.pipe_controler import *
from eModules.fasta_operator import *
import logging
import shutil
from Bio import SeqIO
import collections

usage = """python eProbe_SNP_preprocessor.py -v *.vcf.gz -r ref.fasta -t threads -o output"""
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-v", "--vcf", dest="vcf", action="store", nargs='?',
                    help="Gzipped and indexed VCF file (vcf.gz).",
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
parser.add_argument("--indel_filter", dest="indel_filter", action="store", nargs='?',
                    help="Filter out those SNPs around structural variants in the VCF (on/off, default: on).",
                    choices=["on", "off"], default="on", metavar="STRING")
parser.add_argument("--keep_bed", dest="keep_bed", action="store", nargs='?',
                    help="Genomic regions (BED) that need to retain their SNPs.",
                    metavar="FILE")
parser.add_argument("--keep_distance", dest="keep_distance", action="store", nargs='?',
                    help="Distance flanking(+)/within(-) the regions needed to be kept (default: -50). "
                         "The",
                    default="-50", metavar="STRING")
parser.add_argument("--remove_bed", dest="remove_bed", action="store", nargs='?',
                    help="Genome regions that need to retain their SNPs (bed format)."
                         "In the scenario where both keep_bed and remove_bed are provided as inputs, "
                         "the program will first retain the SNPs in keep_bed, "
                         "and subsequently filter out the remaining SNPs that are located within remove_bed.",
                    metavar="FILE")
parser.add_argument("--rm_distance", action="store", nargs='?', default=50,
                    help="Distance around(+)/within(-) the regions needed to be removed (default: 50).",
                    metavar="STRING")
parser.add_argument("--cluster_flank", dest="cluster_flank", action="store", nargs='?',
                    help="The flank length for filtering out clustered SNPs (default: 100bp).",
                    default="100", metavar="STRING")
parser.add_argument("--max_cluster_snp", dest="max_cluster_snp", action="store", nargs='?',
                    help="Maximum number in flanking regions (default: 3).",
                    default="3", metavar="STRING")
parser.add_argument("--min_cluster_snp", dest="min_cluster_snp", action="store", nargs='?',
                    help="Minimum number in flanking regions (default: 1).",
                    default="1", metavar="STRING")

args = parser.parse_args()

# required params
OUTPUT_PREFIX = args.output
VCF = args.vcf
REFERENCE = args.reference
THREAD = args.thread

# optional params
INDEL_FILTER = args.indel_filter
KEEP_BED = args.keep_bed
REMOVE_BED = args.remove_bed
KEEP_DISTANCE = args.keep_distance
REMOVE_DISTANCE = args.remove_distance
CLUSTER_FLANK = args.cluster_flank
MAX_CLUSTER_SNP = args.max_cluster_snp
MIN_CLUSTER_SNP = args.min_cluster_snp

# make script searching list
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "eModules")
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

# make log file
logging.basicConfig(filename=f"{OUTPUT_PREFIX}.preprocessor.log", level=logging.INFO)

# generate chr_len.tsv
generate_chr_length_tsv(read_fasta(REFERENCE), OUTPUT_PREFIX)

# run Get_variants_from_VCF.py

get_variants_script_name = "Get_variants_from_VCF.py"
get_variants_params = ["-g", f"{OUTPUT_PREFIX}_chr_len.tsv", "-v", VCF, "--indel_filter", INDEL_FILTER, "-t",
                       THREAD, "-o", OUTPUT_PREFIX]
run_python_script(get_variants_script_name, PATH_SEARCH_LIST, get_variants_params)

# run BED_SNPs_filter.py
if REMOVE_BED or KEEP_BED:
    if INDEL_FILTER == "on" and REMOVE_BED:
        merge_command = f"cat {REMOVE_BED} {OUTPUT_PREFIX}.processing_InDels.tsv > {REMOVE_BED}.temp"
        subprocess.run(merge_command, shell=True)
        REMOVE_BED = f"{REMOVE_BED}.temp"
    elif INDEL_FILTER == "on" and not REMOVE_BED:
        REMOVE_BED = f"{OUTPUT_PREFIX}.processing_InDels.tsv"
    else:
        REMOVE_BED = REMOVE_BED

    os.renames(f"{OUTPUT_PREFIX}.processing_SNPs.tsv", f"{OUTPUT_PREFIX}.temp_SNPs.tsv")
    bed_filter_script_name = "BED_SNPs_filter.py"
    bed_filter_params = ["-s", f"{OUTPUT_PREFIX}.temp_SNPs.tsv", "-g", f"{OUTPUT_PREFIX}_chr_len.tsv", "-k",
                         KEEP_BED, "-d", KEEP_DISTANCE, "-r", REMOVE_BED, "-m", REMOVE_DISTANCE, "-o", OUTPUT_PREFIX]
    run_python_script(bed_filter_script_name, PATH_SEARCH_LIST, bed_filter_params)
    os.remove(f"{OUTPUT_PREFIX}.temp_SNPs.tsv")
    if os.path.exists(f"{REMOVE_BED}.temp"):
        os.remove(f"{REMOVE_BED}.temp")

# 3 run Clustered_SNPs_filter.py
cluster_filter_script_name = "Clustered_SNPs_filter.py"
os.renames(f"{OUTPUT_PREFIX}.processing_SNPs.tsv", f"{OUTPUT_PREFIX}.temp_SNPs.tsv")
cluster_filter_params = ["-s", f"{OUTPUT_PREFIX}.temp_SNPs.tsv", "-f", CLUSTER_FLANK, "--max",
                         MAX_CLUSTER_SNP, "-min", MIN_CLUSTER_SNP, "-t", THREAD, "-o", OUTPUT_PREFIX]
run_python_script(cluster_filter_script_name, PATH_SEARCH_LIST, cluster_filter_params)
os.remove(f"{OUTPUT_PREFIX}.temp_SNPs.tsv")
os.rename(f"{OUTPUT_PREFIX}.processing_SNPs.tsv", f"{OUTPUT_PREFIX}.preprocessed_SNPs.tsv")
