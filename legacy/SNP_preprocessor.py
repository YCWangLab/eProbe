import argparse
import os
import logging
import fileinput

from eModules.pipe_controler import *
from eModules.fasta_operator import *

usage = '''python SNP_preprocessor.py -v vcf.gz -r ref.fasta -o output -t threads '''
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-v', '--vcf', dest='vcf', action='store', nargs='?',
                    help='Gzipped and indexed VCF file (vcf.gz).',
                    metavar='FILE', required=True)
parser.add_argument('-r', '--reference', dest='reference', action='store', nargs='?',
                    help='The reference genome (fasta) used to generate probes. '
                         'The same one used for generating VCF.',
                    metavar='FILE', required=True)
parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                    help='Output prefix for temp and result files.',
                    metavar='STRING', required=True)
parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                    help='Number of thread. The number of chromosomes/contigs is recommended (default: 1).'
                         'Multiprocessing will only be carried out on a per-chromosome/contig basis.',
                    default='1', metavar='STRING')
parser.add_argument('--keep_bed', dest='keep_bed', action='store', nargs='?',
                    help='Genomic regions (BED) that need to retain their SNPs.',
                    metavar='FILE')
parser.add_argument('--keep_distance', dest='keep_distance', action='store', nargs='?',
                    help='Distance overhang(+)/recess(-) the regions needed to be kept (default: -50). ',
                    default='-50', metavar='STRING')
parser.add_argument('--remove_bed', dest='remove_bed', action='store', nargs='?',
                    help='Genome regions that need to retain their SNPs (bed format).'
                         'In the scenario where both keep_bed and remove_bed are provided as inputs, '
                         'the program will first retain the SNPs in keep_bed, '
                         'and subsequently filter out the remaining SNPs that are located within remove_bed.',
                    metavar='FILE')
parser.add_argument('--rm_distance', dest='remove_distance', action='store', nargs='?',
                    help='Distance overhang(+)/recess(-) the regions needed to be removed (default: 50).',
                    default='50', metavar='STRING')
parser.add_argument('--cluster_filter', dest='cluster_filter', action='store', nargs='?',
                    help='Enable cluster filter to filter out SNP clusters (on/off, default: on).',
                    choices=['on', 'off'], default='on', metavar='STRING')
parser.add_argument('--cluster_flank', dest='cluster_flank', action='store', nargs='?',
                    help='The flank length for filtering out SNP clusters (default: 60bp).',
                    default='60', metavar='STRING')
parser.add_argument('--max_cluster_snp', dest='max_cluster_snp', action='store', nargs='?',
                    help='Maximum number of SNPs in flanking regions (default: 3).',
                    default='3', metavar='STRING')
parser.add_argument('--min_cluster_snp', dest='min_cluster_snp', action='store', nargs='?',
                    help='Minimum number of SNPs in flanking regions (default: 1).',
                    default='1', metavar='STRING')

args = parser.parse_args()

VCF = args.vcf
REFERENCE = args.reference
OUTPUT_BASE = args.output
THREAD = args.thread

KEEP_BED = args.keep_bed
REMOVE_BED = args.remove_bed
KEEP_DISTANCE = args.keep_distance
REMOVE_DISTANCE = args.remove_distance

CLUMP_FILTER = args.cluster_filter
CLUMP_FLANK = args.cluster_flank
MAX_CLUMP_SNP = args.max_cluster_snp
MIN_CLUMP_SNP = args.min_cluster_snp

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, 'eModules')
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

log_filename = f'{OUTPUT_BASE}.preprocessor.log'
if os.path.exists(log_filename):
    os.remove(log_filename)
logging.basicConfig(filename=log_filename, format='%(message)s', level=logging.INFO)

REFERENCE_DICT = read_fasta(REFERENCE)

if int(THREAD) > len(REFERENCE_DICT.keys()):
    THREAD = str(len(REFERENCE_DICT.keys()))

CHR_SIZE_TSV = get_chr_size(REFERENCE_DICT, OUTPUT_BASE)
TEMP_SNP = f'{OUTPUT_BASE}.temp_SNPs.tsv'
PROCESSING_SNP = f'{OUTPUT_BASE}.processing_SNPs.tsv'

get_snps_script = 'Get_snp_from_vcf.py'
get_snps_params = ['-g', CHR_SIZE_TSV, '-v', VCF, '-t',
                   THREAD, '-o', OUTPUT_BASE]
run_python_script(get_snps_script, PATH_SEARCH_LIST, get_snps_params)

check_file_exists(PROCESSING_SNP, success=f'Successfully retrieved SNPs from {VCF}.',
                  fail=f'Failed to retrieved SNPs from {VCF}.')

if REMOVE_BED or KEEP_BED:
    os.renames(PROCESSING_SNP, TEMP_SNP)
    bed_filter_script = 'Snp_filter_bed.py'
    bed_filter_params = ['-s', TEMP_SNP, '-g', CHR_SIZE_TSV, '-o', OUTPUT_BASE]
    bed_filter_opt_params = {
        '-k': KEEP_BED,
        '-d': KEEP_DISTANCE,
        '-r': REMOVE_BED,
        '-m': REMOVE_DISTANCE
    }
    for param, value in bed_filter_opt_params.items():
        if value:
            bed_filter_params.extend([param, value])
    run_python_script(bed_filter_script, PATH_SEARCH_LIST, bed_filter_params)

    check_file_exists(PROCESSING_SNP, success=f'Successfully filtered SNPs from input BED(s).',
                      fail=f'Failed to filtered SNPs from input BED(s).')

if CLUMP_FILTER == 'on':
    cluster_filter_script = 'Snp_filter_cluster.py'
    os.renames(PROCESSING_SNP, TEMP_SNP)
    cluster_filter_params = ['-s', TEMP_SNP, '-f', CLUMP_FLANK, '--max',
                           MAX_CLUMP_SNP, '--min', MIN_CLUMP_SNP, '-t', THREAD, '-o', OUTPUT_BASE]
    run_python_script(cluster_filter_script, PATH_SEARCH_LIST, cluster_filter_params)

    check_file_exists(PROCESSING_SNP, success=f'Successfully filtered SNP clusters.',
                      fail=f'Failed to filter SNP clusters.')

else:
    first_line = '\t'.join(['chr', 'pos', 'type', 'ref', 'alt']) + '\n'
    with fileinput.FileInput(PROCESSING_SNP, inplace=True) as file:
        for line in file:
            if file.isfirstline():
                print(first_line, end='')
            print(line, end='')

if os.path.exists(TEMP_SNP):
    os.remove(TEMP_SNP)
if os.path.exists(CHR_SIZE_TSV):
    os.remove(CHR_SIZE_TSV)
os.rename(PROCESSING_SNP, f'{OUTPUT_BASE}.preprocessed_SNPs.tsv')
