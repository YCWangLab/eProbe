import argparse
import os
import logging

from eModules.pipe_controler import *
from eModules.fasta_operator import *

usage = '''python SNP_subsampler.py -f snp.tsv -r ref.fasta -t threads -o output'''
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-f', '--snp_df', dest='snp_df', action='store', nargs='?',
                    help='Input SNPs data frame (tsv).',
                    metavar='FILE', required=True)
parser.add_argument('-r', '--reference', dest='reference', action='store', nargs='?',
                    help='The reference genome (fasta) used to generate probes. '
                         'The same one used for VCF calling.',
                    metavar='FILE', required=True)
parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                    help='Number of thread (default: 1).',
                    default='1', metavar='STRING')
parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                    help='Output prefix for process files and result files.',
                    metavar='STRING', required=True)
parser.add_argument('--window_size', dest='window_size', action='store', nargs='?',
                    help='Window size for sub-sampling (default: 1000). '
                         'eProbe will select one SNP from each non-overlapping window containing SNPs.',
                    default='1000', metavar='STRING')
parser.add_argument('--output_mode', dest='output_mode', action='store', nargs='?',
                    help='Output mode (both/feature/window, default: window).'
                         'window: If SNPs in a window overlapping the BED (if provided), prioritizing the one in'
                         'the interval specified in the BED. Otherwise, then either select one based on '
                         'weights (if provided), or randomly choose one.'
                         'feature: Only output all SNPs falling within the BED.'
                         'both: Output both Window and Feature SNP sets and their merged SNP set',
                    choices=['both', 'feature', 'window'], default='window', metavar='STRING')
parser.add_argument('--desire_tm', dest='desire_tm', action='store', nargs='?',
                    help='Desire melting temperature (default: 75).',
                    default='75', metavar='STRING')
parser.add_argument('--desire_gc', dest='desire_gc', action='store', nargs='?',
                    help='Desire GC content (default: 50).',
                    default='50', metavar='STRING')
parser.add_argument('--select_weights', dest='weights', action='store', nargs='?',
                    help='Weights for gc, tm, dust, hairpin, dimer in order.'
                         'eProbe will calculate the score by weighted summing the normalized distance of each'
                         'biophysical/biochemical property to their desired value (for hairpin, dimer and DUST score, '
                         'the smaller the better), then select the SNP with the best overall performance.'
                         'If not provided, eProbe will randomly select one SNP from a window.',
                    default=None, metavar='STRING')
parser.add_argument('--keep_bed', dest='keep_bed', action='store', nargs='?',
                    help='Annotation file. Usually a 4-column BED file converted from gff, '
                         'with last column describing the function of this genomic region. '
                         'If provided, a SNP of a window will be kept prioritarily if it overlaps with BED.',
                    default=None, metavar='FILE')
parser.add_argument('--keep_bed_column', dest='keep_bed_column', action='store', nargs='?',
                    help='Number of the BED file column used to match feature. '
                         'If not provided, keep all genomic regions of the BED file.',
                    default=None, metavar='STRING')
parser.add_argument('--keep_bed_feature', dest='keep_bed_feature', action='store', nargs='?',
                    help='What feature in the column want to be kept (e.g., exon). '
                         'If not provided, keep all genomic regions of the BED file.',
                    default=None, metavar='STRING')
parser.add_argument('--probe_number', dest='probe_number', action='store', nargs='?',
                    help='Randomly subsample SNPs to a certain number.'
                         'If the selected SNPs number is less than the input SNPs number, output selected SNPs.'
                         'If output mode is Both, subsample based on merged SNPs dataframe.',
                    default=None, metavar='STRING')
parser.add_argument('--seed', dest='seed', action='store', nargs='?',
                    help='Seed number for random subsampling.',
                    default='123', metavar='STRING')

args = parser.parse_args()

SNP_DF = args.snp_df
REFERENCE = args.reference
OUTPUT_BASE = args.output
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

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, 'eModules')
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

log_filename = f'{OUTPUT_BASE}.subsampler.log'
if os.path.exists(log_filename):
    os.remove(log_filename)
logging.basicConfig(filename=log_filename, format='%(message)s', level=logging.INFO)

CHR_SIZE_TSV = get_chr_size(read_fasta(REFERENCE), OUTPUT_BASE)

selector_script_name = 'Snp_selector.py'
selector_params = ['-d', SNP_DF, '-c', CHR_SIZE_TSV, '-s', WINDOW_SIZE, '-m', SELECT_MODE, '-r', THREAD, '-o',
                   OUTPUT_BASE]
selector_opt_params = {
    '-t': DESIRE_TM,
    '-g': DESIRE_GC,
    '-w': SELECT_WEIGHTS,
    '-b': KEEP_BED,
    '-n': KEEP_BED_COLUMN,
    '-f': KEEP_BED_FEATURE,
}

for param, value in selector_opt_params.items():
    if value:
        selector_params.extend([param, value])
run_python_script(selector_script_name, PATH_SEARCH_LIST, selector_params)

if SELECT_MODE == 'window':
    SELECTED_SNP = f'{OUTPUT_BASE}_window_SNPs.tsv'
elif SELECT_MODE == 'feature':
    SELECTED_SNP = f'{OUTPUT_BASE}_feature_SNPs.tsv'
elif SELECT_MODE == 'both':
    SELECTED_SNP = f'{OUTPUT_BASE}_merged_SNPs.tsv'
check_file_exists(SELECTED_SNP, success=f'Successfully subsampled SNPs based on the sliding window method.',
                  fail=f'Failed to subsample SNPs based on the sliding window method.')

if PROBE_NUMBER:
    subsampler_script_name = 'Snp_subsampler.py'
    subsampler_params = ['-s', SEED, '-r', PROBE_NUMBER, '--df1', SELECTED_SNP, '-o',
                         OUTPUT_BASE]
    run_python_script(subsampler_script_name, PATH_SEARCH_LIST, subsampler_params)

    check_file_exists(f'{OUTPUT_BASE}.subsampled_{PROBE_NUMBER}.tsv',
                      success=f'Randomly selected {PROBE_NUMBER} SNPs.',
                      fail=f'Failed to select SNPs to a certain size randomly.')
else:
    print('Disable random selection to a certain number of SNPs. Skipping...')

if os.path.exists(CHR_SIZE_TSV):
    os.remove(CHR_SIZE_TSV)

