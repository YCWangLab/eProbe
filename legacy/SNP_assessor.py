import argparse
import os
import logging
import shutil

from eModules.pipe_controler import *
from eModules.fasta_operator import *

usage = '''python SNP_assessor.py -f SNPs_dataframe.tsv --evaluator [options] '''
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-f', '--snp_df', dest='snp_df', action='store', nargs='?',
                    help='SNPs data frame.',
                    metavar='FILE', required=True)
parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                    help='[TAG] Number of thread (default: 1).',
                    default='1', metavar='STRING')
parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                    help='Output prefix.', metavar='STRING', required=True)
parser.add_argument('-a', '--assessor', dest='assessor', action='store', nargs='?',
                    help='By what do you want to assess the probe set (tag/distance/missing).'
                         'Tag: plotting the distribution of different biophysical/biochemical tags '
                         'to evaluate the performance in hybridization capture.'
                         'Distance: calculate pairwise distance matrix based on the SNPs covered by the probe'
                         'set and compare it with that of the original genome-wide data.'
                         'Missing: simulate random loss of the probe-covered SNPs and perform '
                         'PCA to evaluate the tolerance to missing rate of genotyping.',
                    choices=['tag', 'distance', 'missing'], metavar='FILE', required=True)
parser.add_argument('--tag', dest='tag', action='store', nargs='?',
                    help='[TAG] Which biophysical/biochemical feautre(s) [gc/tm/complexity/hairpin/dimer, separated by comma] is/are used for '
                         'drawing distribution plots (default: all)?',
                    default='gc,tm,complexity,hairpin,dimer', metavar='STRING')
parser.add_argument('-r', '--reference', dest='reference', action='store', nargs='?',
                    help='[TAG] The reference genome (fasta) used to generate probes.'
                         'The same one used for VCF calling.',
                    metavar='STRING')
parser.add_argument('-v', '--vcf', dest='vcf', action='store', nargs='?',
                    help='[DISTANCE and MISSING] Gzipped and indexed VCF file (vcf.gz).'
                         'Required for evaluating with distance and missing.',
                    metavar='FILE')
parser.add_argument('--stat', dest='TAG_stat', action='store', nargs='?',
                    help='[TAG] Stat style [percent/density/probability/count/frequency] for plot (default: percent).',
                    choices=['percent', 'density', 'probability', 'count', 'frequency'], default='percent',
                    metavar='STRING')
parser.add_argument('--bins', dest='TAG_bins', action='store', nargs='?',
                    help='[TAG] Bin interval for histogram (default: auto).',
                    default='auto', metavar='STRING')
parser.add_argument('--xlim', dest='TAG_xlim', action='store', nargs='?',
                    help='[TAG] Lower and upper limits of X axis, separated by comma (default: auto).',
                    default='auto', metavar='STRING')
parser.add_argument('--ylim', dest='TAG_ylim', action='store', nargs='?',
                    help='[TAG] Lower and upper limits of Y axis, separated by comma (default: auto).',
                    default='auto', metavar='STRING')
parser.add_argument('--colors', dest='colors', action='store', nargs='?',
                    help='[TAG] Colors for plot distribution (default: 5 default colors).',
                    default=None, metavar='STRING')
parser.add_argument('--correlation', dest='correlation', action='store', nargs='?',
                    help='[DISTANCE] Method for calculating the correlation coefficient (default: Pearson).',
                    choices=['Pearson', 'Spearman'], default='Pearson', metavar='STRING')
parser.add_argument('--cluster', dest='cluster', action='store', nargs='?',
                    help='[DISTANCE] Conduct hierarchical clustering (on/off, default: on).',
                    default='on', metavar='STRING')
parser.add_argument('--cmap', dest='cmap', action='store', nargs='?',
                    help='[DISTANCE] Cmap for plotting heatmap (try flare/crest/viridis_r, default: YlOrBr). ',
                    default='YlOrBr', metavar='STRING')
parser.add_argument('--missing_rate', dest='missing_rate', action='store', nargs='?',
                    help='[MISSING] Missing rate. Can be fixed for all sites (e.g., 0.2), '
                         'or a list for random selection (default: 0.2,0.4,0.6)',
                    default='0.2,0.4,0.6', metavar='STRING')
parser.add_argument('--missing_genotype', dest='missing_genotype', action='store', nargs='?',
                    help='[MISSING] Code for missing genotype, normally \'./.:0,0,0\' for VCF called from BCFtools'
                         'and \'.:0,0:0:.:0,0\' for GATK)',
                    default='./.:0,0,0', metavar='STRING')
parser.add_argument('--pc', dest='pcn', action='store', nargs='?',
                    help='[MISSING] How many PCs to draw (default: 1,2).',
                    default='1,2', metavar='FILE')
parser.add_argument('--pop_data', dest='pop_data', action='store', nargs='?',
                    help='[MISSING] A population information to color the map.'
                         'Should be a tab-separated file and each line consist of '
                         'Sample_ID and Population_ID, without header.',
                    default=None, metavar='FILE')
parser.add_argument('--TAG_tm_method', dest='TAG_tm_method', action='store', nargs='?',
                    help='[TAG] Method of calculating melting temperature (tm_NN/tm_GC; default: tm_NN).',
                    choices=['tm_NN', 'tm_GC'], default='tm_NN', metavar='STRING')
parser.add_argument('--TAG_tm_table', dest='TAG_tm_table', action='store', nargs='?',
                    help='[TAG] Table of thermodynamic NN values.'
                         'DNA/DNA hybridizations:'
                         'DNA_NN1: values from Breslauer et al. (1986)'
                         'DNA_NN2: values from Sugimoto et al. (1996)'
                         'DNA_NN3: values from Allawi & SantaLucia (1997)'
                         'DNA_NN4: values from SantaLucia & Hicks (2004)'
                         'For RNA/DNA hybridizations:'
                         'R_DNA_NN1: values from Sugimoto et al. (1995) (default).',
                    choices=['DNA_NN1', 'DNA_NN2', 'DNA_NN3', 'DNA_NN4', 'R_DNA_NN1'],
                    default='R_DNA_NN1', metavar='STRING')
parser.add_argument('--TAG_hp_ms', dest='TAG_match_score', action='store', nargs='?',
                    help='[TAG] Match score for calculating hairpin score (default: 1).',
                    default='1', metavar='STRING')
parser.add_argument('--TAG_hp_mmp', dest='TAG_mismatch_penalty', action='store', nargs='?',
                    help='[TAG] Mismatch penalty for calculating hairpin score (default: -1).',
                    default='-1', metavar='STRING')
parser.add_argument('--TAG_hp_gop', dest='TAG_gapopen_penalty', action='store', nargs='?',
                    help='[TAG] Gap open penalty for calculating hairpin score (default: -5).',
                    default='-5', metavar='STRING')
parser.add_argument('--TAG_hp_gep', dest='TAG_gapextend_penalty', action='store', nargs='?',
                    help='[TAG] Gap extension penalty for calculating hairpin score (default: -5).',
                    default='-5', metavar='STRING')
parser.add_argument('--TAG_dm_kmer', dest='TAG_kmer_length', action='store', nargs='?',
                    help='[TAG] Length of kmer for calculating dimer score (default: 11).',
                    default='11', metavar='STRING')
parser.add_argument('--TAG_dm_min_freq', dest='TAG_kmer_minfreq', action='store', nargs='?',
                    help='[TAG] Minimum value of kmer frequency for calculating dimer score (default: 2).',
                    default='2', metavar='STRING')
args = parser.parse_args()

SNP_DF = args.snp_df
REFERENCE = args.reference
VCF = args.vcf
THREAD = args.thread
OUTPUT_BASE = args.output
EVALUATOR = args.assessor
COLORS = args.colors

if args.tag == 'all':
    TAG = 'gc,tm,complexity,hairpin,dimer'
else:
    TAG = args.tag
TAG_STAT = args.TAG_stat
TAG_BINS = args.TAG_bins
TAG_XLIM = args.TAG_xlim
TAG_YLIM = args.TAG_ylim
TAG_TM_METHOD = args.TAG_tm_method
TAG_TM_TABLE = args.TAG_tm_table
TAG_HAIRPIN_MATCH_SCORE = args.TAG_match_score
TAG_HAIRPIN_MISMATCH_PENALTY = args.TAG_mismatch_penalty
TAG_HAIRPIN_GAPOPEN_PENALTY = args.TAG_gapopen_penalty
TAG_HAIRPIN_GAPEXTEND_PENALTY = args.TAG_gapextend_penalty
TAG_DIMER_KMER = args.TAG_kmer_length
TAG_DIMER_MIN_FREQ = args.TAG_kmer_minfreq

CLUSTER = args.cluster
CMAP = args.cmap
CORRELATION = args.correlation

MISSING_RATE = args.missing_rate
MISSING_GENOTYPE = args.missing_genotype
PCN = args.pcn
POP_DATA = args.pop_data

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, 'eModules')
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

log_filename = f'{OUTPUT_BASE}.evaluator.log'
if os.path.exists(log_filename):
    os.remove(log_filename)
logging.basicConfig(filename=log_filename, format='%(message)s', level=logging.INFO)

if EVALUATOR.lower() == 'tag':
    print('Start assessing biophysical/biochemical features...')

    get_probes_script = 'Get_probe_from_snp.py'
    get_probes_params = ['-d', SNP_DF, '-g', REFERENCE, '-m', 'edge', '-o', OUTPUT_BASE]
    run_python_script(get_probes_script, PATH_SEARCH_LIST, get_probes_params)

    PROCESSING_PROBE = f'{OUTPUT_BASE}.filtered_probe.fasta'

    get_tag_gc_script = 'Get_tag_gc.py'
    get_tag_tm_script = 'Get_tag_tm.py'
    get_tag_dust_script = 'Get_tag_complexity.py'
    get_tag_hairpin_script = 'Get_tag_hairpin.py'
    get_tag_dimer_script = 'Get_tag_dimer.py'

    tmp_dir = f'{OUTPUT_BASE}_eProbe_temp'
    make_temp_dir(tmp_dir)

    tmp_output_base = os.path.join(tmp_dir, os.path.basename(OUTPUT_BASE))

    if 'gc' in TAG.split(','):
        gc_script_params = ['-f', PROCESSING_PROBE, '-t', THREAD, '-o',
                            f'{tmp_output_base}.gc']
        run_python_script(get_tag_gc_script, PATH_SEARCH_LIST, gc_script_params)

    if 'tm' in TAG.split(','):
        tm_script_params = ['-f', PROCESSING_PROBE, '-t', THREAD, '-m', TAG_TM_METHOD, '--table', TAG_TM_TABLE,
                            '-o', f'{tmp_output_base}.tm']
        run_python_script(get_tag_tm_script, PATH_SEARCH_LIST, tm_script_params)

    if 'complexity' in TAG.split(','):
        complexity_script_params = ['-f', PROCESSING_PROBE, '-t', THREAD, '-o',
                                    f'{tmp_output_base}.complexity']
        run_python_script(get_tag_dust_script, PATH_SEARCH_LIST, complexity_script_params)

    if 'hairpin' in TAG.split(','):
        hairpin_script_params = ['-f', PROCESSING_PROBE, '-t', THREAD,
                                 '--mp', TAG_HAIRPIN_MATCH_SCORE, '--mmp', TAG_HAIRPIN_MISMATCH_PENALTY, '--gp',
                                 TAG_HAIRPIN_GAPOPEN_PENALTY, '--gpe', TAG_HAIRPIN_GAPEXTEND_PENALTY,
                                 '-o', f'{tmp_output_base}.hairpin']
        run_python_script(get_tag_hairpin_script, PATH_SEARCH_LIST, hairpin_script_params)

    if 'dimer' in TAG.split(','):
        dimer_script_params = ['-f', PROCESSING_PROBE, '-k', TAG_DIMER_KMER, '--min_freq',
                               TAG_DIMER_MIN_FREQ, '-t', THREAD, '-o',
                               f'{tmp_output_base}.dimer']
        run_python_script(get_tag_dimer_script, PATH_SEARCH_LIST, dimer_script_params)

    selected_tags = get_unique_file_extensions(tmp_dir)
    if len(selected_tags):
        merge_tags_script = 'Tags_merger.py'
        merge_tags_script_params = ['-o', tmp_output_base, '-f']

        for tag in selected_tags:
            if tag in ['.tm', '.gc', '.complexity', '.hairpin', '.dimer']:
                merge_tags_script_params = merge_tags_script_params + [f'{tmp_output_base}{tag}']
        run_python_script(merge_tags_script, PATH_SEARCH_LIST, merge_tags_script_params)

    draw_script = 'Draw_dist_from_df.py'
    draw_params = ['-d', f'{tmp_output_base}.merged_tags.tsv', '--stat', TAG_STAT, '--bins', TAG_BINS, '--xlim',
                   TAG_XLIM,
                   '--ylim', TAG_YLIM, '-c', ' '.join(TAG.split(',')), '-o', OUTPUT_BASE]
    if COLORS:
        draw_params = draw_params + ['--colors', COLORS]
    run_python_script(draw_script, PATH_SEARCH_LIST, draw_params)

    output_files = []
    for tag in selected_tags:
        if tag in ['.tm', '.gc', '.complexity', '.hairpin', '.dimer']:
            output_files.append(f'{OUTPUT_BASE}_{tag.split(".")[1]}_dist.jpg')
    check_file_exists(output_files, success=f'Accessing biophysical/biochemical features completed.',
                      fail=f'Accessing biophysical/biochemical features failed.')

    shutil.rmtree(tmp_dir)
    os.remove(PROCESSING_PROBE)

elif EVALUATOR.lower() == 'distance':
    get_vcf_script = 'Get_vcf_from_snp.py'
    get_vcf_params = ['-f', SNP_DF, '-v', VCF, '-o', OUTPUT_BASE]
    run_python_script(get_vcf_script, PATH_SEARCH_LIST, get_vcf_params)

    probe_vcf = f'{OUTPUT_BASE}.vcf'

    get_distance_script = 'Assess_by_distance.py'
    get_distance_params = ['--vcf1', VCF, '--vcf2', probe_vcf, '--cluster', CLUSTER, '--cor', CORRELATION, '--cmap',
                           CMAP, '-o',
                           OUTPUT_BASE]
    run_python_script(get_distance_script, PATH_SEARCH_LIST, get_distance_params)

    output_files = [f'{OUTPUT_BASE}_All_distance.jpg', f'{OUTPUT_BASE}_Part_distance.jpg']
    check_file_exists(output_files, success=f'Accessing estimated pairwise distance completed.',
                      fail=f'Accessing estimated pairwise distance failed.')


elif EVALUATOR.lower() == 'missing':
    get_vcf_script = 'Get_vcf_from_snp.py'
    get_vcf_params = ['-f', SNP_DF, '-v', VCF, '-o', OUTPUT_BASE]
    run_python_script(get_vcf_script, PATH_SEARCH_LIST, get_vcf_params)

    probe_vcf = f'{OUTPUT_BASE}.vcf'

    get_missing_script = 'Assess_by_missing.py'
    get_missing_params = ['-v', probe_vcf, '-r', MISSING_RATE, '-g', MISSING_GENOTYPE, '--pc', PCN, '-o', OUTPUT_BASE]
    if POP_DATA:
        get_missing_params = get_missing_params + ['--pop_data', POP_DATA]
    run_python_script(get_missing_script, PATH_SEARCH_LIST, get_missing_params)

    import itertools
    PCs = list(itertools.combinations(PCN.split(','), 2))
    output_files = []
    for PC in PCs:
        output_files.append(f'''{OUTPUT_BASE}_original_{f'PC{PC[0]}'}_{f'PC{PC[1]}'}.jpg''')
    check_file_exists(output_files, success=f'Accessing estimated pairwise distance completed.',
                      fail=f'Accessing estimated pairwise distance failed.')

else:
    raise ValueError('Choose your evaluator: tag, distance, missing')
