import argparse
import logging
import os

from eModules.pipe_controler import *
from eModules.fasta_operator import *

usage = '''python Seq_generator.py --fasta *.fasta --length 100 --step 30 -o output'''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('-m', '--method', dest='method', action='store', nargs='?',
                    help='Method for generating probes with sequences (seq/bed, default: seq). '
                         'seq: generating probes with fasta file. '
                         'bed: generating probes with BED file and corresponding reference genome.',
                    choices=['seq', 'bed'], default='seq', metavar='STRING')
parser.add_argument('--haplotyping', dest='haplotyping', action='store', nargs='?',
                    help='Infer haplotype sequences of each interval specified by BED file using VCF file'
                         'or input sequences of different alleles indicated (on/off, default: off).',
                    choices=['on', 'off'], default='off', metavar='STRING')
parser.add_argument('-l', '--length', dest='length', action='store', nargs='?',
                    help='Length of probes.',
                    metavar='STRING', required=True)
parser.add_argument('-s', '--step', dest='step', action='store', nargs='?',
                    help='Step for generating probes using sliding windows.',
                    metavar='STRING', required=True)
parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                    help='Number of thread (default: 1).',
                    default='1', metavar='STRING')
parser.add_argument('-f', '--fasta', dest='fasta', action='store', nargs='?',
                    help='[SEQ] Input template sequences (*.fasta) for generating probes.',
                    metavar='FILE')
parser.add_argument('-b', '--bed', dest='bed', action='store', nargs='?',
                    help='[BED] Genomic regions for subtracting sequences and haplotyping (BED format).',
                    metavar='FILE')
parser.add_argument('-r', '--reference', dest='reference', action='store', nargs='?',
                    help='[BED] The reference genome used to generate probes with BED file. ',
                    metavar='FILE')
parser.add_argument('-v', '--vcf', dest='vcf', action='store', nargs='?',
                    help='[BED] Input vcf.gz (with .tbi) for phasing with shapeit (should be in your $PATH).',
                    metavar='FILE')
parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                    help='Output prefix.',
                    metavar='STRING', required=True)
parser.add_argument('--sep', dest='sep', action='store', nargs='?',
                    help='[SEQ] Separator for indicating allele num (default: "_"). ',
                    default='_', metavar='STRING')
parser.add_argument('--freq', dest='freq', action='store', nargs='?',
                    help='[BED] Minimum frequency of retained alleles (default: 0.05) in haplotyping.',
                    default='0.05', metavar='STRING')


args = parser.parse_args()

METHOD = args.method
HAPLOTYPING = args.haplotyping
PROBE_LENGTH = args.length
PROBE_STEP = args.step
THREAD = args.thread

FASTA = args.fasta
SEPARATOR = args.sep

BED = args.bed
REFERENCE = args.reference
VCF = args.vcf
MINIMUM_FREQUENCY = args.freq

OUTPUT_BASE = args.output

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, 'eModules')
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

log_filename = f'{OUTPUT_BASE}.Seq_generator.log'
if os.path.exists(log_filename):
    os.remove(log_filename)
logging.basicConfig(filename=log_filename, format='%(message)s', level=logging.INFO)

if METHOD.lower() == 'seq':
    if HAPLOTYPING == 'on':
        generator_script_name = 'Get_probe_from_allele.py'
        generator_params = ['--input', METHOD, '-f', FASTA, '-l', PROBE_LENGTH, '-s', PROBE_STEP, '--min_freq',
                            MINIMUM_FREQUENCY, '--sep', SEPARATOR, '-o', OUTPUT_BASE]
        run_python_script(generator_script_name, PATH_SEARCH_LIST, generator_params)

    else:
        generator_script_name = 'Get_probe_from_seq.py'
        generator_params = ['-f', FASTA, '-l', PROBE_LENGTH, '-s', PROBE_STEP, '-o', OUTPUT_BASE]
        run_python_script(generator_script_name, PATH_SEARCH_LIST, generator_params)


elif METHOD.lower() == 'bed':
    if HAPLOTYPING == 'on':
        get_alleles_script_name = 'Get_allele_from_vcf.py'
        get_alleles_params = ['-b', BED, '-r', REFERENCE, '-v', VCF, '-t', THREAD,
                              '-f', MINIMUM_FREQUENCY, '-o', OUTPUT_BASE]
        run_python_script(get_alleles_script_name, PATH_SEARCH_LIST, get_alleles_params)

        check_file_exists(f'{OUTPUT_BASE}.allele.fasta', success='Successfully generated allele sequences.',
                          fail='Failed to generate allele sequences.')

        generator_script_name = 'Get_probe_from_allele.py'
        generator_params = ['--input', METHOD, '-f', f'{OUTPUT_BASE}.allele.fasta', '-l', PROBE_LENGTH, '-s',
                            PROBE_STEP, '--min_freq',
                            MINIMUM_FREQUENCY, '--sep', SEPARATOR, '-o', OUTPUT_BASE]
        run_python_script(generator_script_name, PATH_SEARCH_LIST, generator_params)
        os.remove(f'{OUTPUT_BASE}.allele.fasta')
    else:
        generator_script_name = 'Get_probe_from_bed.py'
        generator_params = ['-b', BED, '-r', REFERENCE, '-l', PROBE_LENGTH, '-s', PROBE_STEP, '-o', OUTPUT_BASE]
        run_python_script(generator_script_name, PATH_SEARCH_LIST, generator_params)


else:
    raise ValueError(f'Method must be either seq or bed. Unvalid {METHOD} detected.')

check_file_exists(f'{OUTPUT_BASE}.probe.fasta', success='Successfully generated probe sequences.',
                  fail='Failed to generate probe sequences.')
