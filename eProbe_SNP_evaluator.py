import argparse
import os
import logging

usage = """python Get_VCF_from_SNPs.py -f SNPs_dataframe.tsv -v raw.vcf.gz(*.tbi) -o out.vcf """
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-f', '--snp_df', dest="snp_df", action="store", nargs='?',
                    help="SNPs dataframe.", metavar="FILE", required=True)
parser.add_argument("--evaluator", dest="evaluator", action="store", nargs='?',
                    help="By what do we evaluate the performance of the probe (tag/distance/missing)."
                         "Tag: plotting the distribution of different biophysical tags "
                         "to evaluate the performance in capture experiments."
                         "Distance: calculate pairwise distance matrix based on the SNPs covered by the probe"
                         "and compare it with that of the original data."
                         "Missing: Simulate random loss of the probe-covered SNPs and perform "
                         "PCA to evaluate the genotype tolerance to missing rates.",
                    choices=["tag", "distance", "missing"], metavar="FILE", required=True)
parser.add_argument("-v", "--vcf", dest="vcf", action="store", nargs='?',
                    help="Gzipped and indexed VCF file (vcf.gz)."
                         "Required for evaluating with distance and missing.",
                    metavar="FILE")

parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output prefix.", metavar="STRING")
args = parser.parse_args()

SNP_DFS = args.df
VCF = args.vcf
OUTPUT_PREFIX = args.output
EVALUATOR = args.evaluator
# tag
TAG = args.tag
STAT = args.stat
BIN = args.bins
# distance
DISTANCE = args.distance
CLUSTER = args.cluster
CMAP = args.cmap
CORRELATION = args.correlation
# missing
MISSING_RATE = args.missing_rate
MISSING_GENOTYPE = args.missing_genotype

# make script searching list
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "eModules")
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

# make log file
logging.basicConfig(filename=f"{OUTPUT_PREFIX}.evaluator.log", level=logging.INFO)

if EVALUATOR.lower() == "tag":
    # 1 read SNPs dataframe
    # 2 extract sequence
    # 3 run tags
    # 4 draw dist
    pass
elif EVALUATOR.lower() == "distance":
    # 1 read SNPs dataframe
    # 2 extract SNPs
    # 3 calculate distance
    # 4 draw
    pass
elif EVALUATOR.lower() == "missing":
    # 1 read SNPs dataframe
    # 2 extract SNPs
    # 3 add missing
    # 4 run pca
    # 5 draw
    pass
else:
    raise ValueError("Choose your evaluator: tag, distance, missing")
