import os
import argparse
import logging
import shutil

from eModules.pipe_controler import *
from eModules.fasta_operator import *

usage = '''python SNP_filter.py -f snp.tsv -r ref.fasta --BG_db BGdb1,BGdb2,BGdb3
 --AC_db ACdb1,ACdb2,ACdb3 --TX_db TXdb1,TXdb2,TXdb3 -t threads -o output'''
parser = argparse.ArgumentParser(description=usage)

parser.add_argument('-f', '--snp_df', dest='snp_df', action='store', nargs='?',
                    help='Input SNPs data frame (tsv).',
                    metavar='FILE', required=True)
parser.add_argument('-r', '--reference', dest='reference', action='store', nargs='?',
                    help='The reference genome (fasta) used to generate probes. '
                         'The same one used for generating VCF.',
                    metavar='FILE', required=True)
parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                    help='Output prefix for process files and result files.',
                    metavar='STRING', required=True)
parser.add_argument('-t', '--thread', dest='thread', action='store', nargs='?',
                    help='Number of thread (default: 1).',
                    default='1', metavar='STRING')
parser.add_argument('-l', '--length', dest='length', action='store', nargs='?',
                    help='Length of probes (default: 81).',
                    default='81', metavar='STRING')
parser.add_argument('-s', '--shift', dest='shift', action='store', nargs='?',
                    help='Position shift (default: 0) fo method pos. - for shifting left '
                         'and + for shifting right.',
                    default='0', metavar='STRING')

parser.add_argument('--BG_db', dest='BGfilter_database', action='store', nargs='?',
                    help='Kraken2 databases (separated with comma) for background noise filtering.',
                    default=None, metavar='STRING')
parser.add_argument('--BG_min_hits', dest='BGfilter_min_hits', action='store', nargs='?',
                    help='Minimum number of non-overlapping kmer hits (separate kmer stacks)'
                         'to filter out (default: 2).',
                    default='2', metavar='STRING')

parser.add_argument('--AC_db', dest='ACfilter_database', action='store', nargs='?',
                    help='Bowtie2 databases (separated with comma) for filtering SNPs in inaccessible regions of genome.',
                    default=None, metavar='FILE')
parser.add_argument('--AC_keep_hit', dest='ACfilter_keep_hit', action='store', nargs='?',
                    help='The number of hits to be kept in bowtie2 mapping in accessibility '
                         'filtering (default: 100).',
                    default='100', metavar='STRING')
parser.add_argument('--AC_map_mode', dest='ACfilter_map_mode', action='store', nargs='?',
                    help='Bowtie2 local mapping mode in accessibility filtering. '
                         '(very-fast/fast/sensitive/very-sensitive, default: sensitive)',
                    default='sensitive', metavar='STRING')
parser.add_argument('--AC_trim5', dest='ACfilter_trim5', action='store', nargs='?',
                    help='Trim bases from 5\'/left end of reads (default: 0).',
                    default='0', metavar='STRING')
parser.add_argument('--AC_trim3', dest='ACfilter_trim3', action='store', nargs='?',
                    help='Trim bases from 3\'/right end of reads (default: 0).',
                    default='0', metavar='STRING')
parser.add_argument('--AC_filter_mode', dest='ACfilter_mode', action='store', nargs='?',
                    help='Accessibility filtering mode (strict/moderate/custom, default: moderate).'
                         'If custom is enabled, custom filter mode parameters will be used.',
                    default='moderate', metavar='STRING')
parser.add_argument('--AC_min_mapq', dest='ACfilter_min_mapq', action='store', nargs='?',
                    help='[custom filter mode parameter] Minimum mapping quality (default: 30).',
                    default='30', metavar='STRING')
parser.add_argument('--AC_max_nm', dest='ACfilter_max_nm', action='store', nargs='?',
                    help='[custom filter mode parameter] Maximum edit distances of alignment (default: 3). ',
                    default='3', metavar='STRING')
parser.add_argument('--AC_max_xo', dest='ACfilter_max_xo', action='store', nargs='?',
                    help='[custom filter mode parameter] Maximum gap number of alignment (default: 0).',
                    default='0', metavar='STRING')
parser.add_argument('--AC_max_xg', dest='ACfilter_max_xg', action='store', nargs='?',
                    help='[custom filter mode parameter] Maximum total gap length of alignment (default: 0).',
                    default='0', metavar='STRING')
parser.add_argument('--AC_max_xsr', dest='ACfilter_max_xsr', action='store', nargs='?',
                    help='[custom filter mode parameter] Maximum secondary alignment score ratio '
                         'compared to best alignment (default: 0.8).',
                    default='0.8', metavar='STRING')

parser.add_argument('--TX_db', dest='TXfilter_database', action='store', nargs='?',
                    help='Bowtie2 databases (separated with comma) for taxonomic filtering.',
                    default=None, metavar='STRING')
parser.add_argument('--TX_keep_hit', dest='TXfilter_keep_hit', action='store', nargs='?',
                    help='The number of hits to be kept in bowtie2 mapping (default: 100).',
                    default='100', metavar='STRING')
parser.add_argument('--TX_taxa', dest='TXfilter_taxa', action='store', nargs='?',
                    help='Taxid of target species/taxa to be extracted from taxonomic identification',
                    metavar='STRING')
parser.add_argument('--TX_names_dmp', dest='TXfilter_names_dmp', action='store', nargs='?',
                    help='NCBI taxonomy names.dmp.', metavar='FILE')
parser.add_argument('--TX_nodes_dmp', dest='TXfilter_nodes_dmp', action='store', nargs='?',
                    help='NCBI taxonomy nodes.dmp.', metavar='FILE')
parser.add_argument('--TX_acc2tax', dest='TXfilter_acc2tax', action='store', nargs='?',
                    help='Accession2taxid.txt, the match file of genome accession and taxid.', metavar='FILE')
parser.add_argument('--TX_minedit', dest='TXfilter_min_edit', action='store', nargs='?',
                    help='Minimum edit distance to assign a sequence in ngsLCA (default: 0).',
                    default='0', metavar='STRING')
parser.add_argument('--TX_maxedit', dest='TXfilter_max_edit', action='store', nargs='?',
                    help='Maximum edit distance to assign a sequence in ngsLCA. All sequences have edit distances '
                         'with reference genome lower than this value are regarded as the same distance (default: 2).',
                    default='2', metavar='STRING')

parser.add_argument('--BFfilter', dest='BFfilter', action='store', nargs='?',
                    help='Enable biophysical filter to filter out probes with undesirable biophysical properties',
                    choices=['on', 'off'], default='off', metavar='STRING')
parser.add_argument('--BF_tm', dest='BFfilter_keep_tm', action='store', nargs='?',
                    help='Upper and lower limits for filtering with melting temperature (default: 68,78).'
                         'Input "off" to disable filtering by melting temperature.',
                    default='60,80', metavar='STRING')
parser.add_argument('--BF_gc', dest='BFfilter_keep_gc', action='store', nargs='?',
                    help='Upper and lower limits for filtering with GC content (default: 40,60).'
                         'Input "off" to disable filtering by GC content.',
                    default='40,60', metavar='STRING')
parser.add_argument('--BF_hairpin', dest='BFfilter_keep_hairpin', action='store', nargs='?',
                    help='Upper and lower limits for filtering with hairpin score (default: 0,60).'
                         'Input "off" to disable filtering by hairpin score.',
                    default='0,60', metavar='STRING')
parser.add_argument('--BF_complexity', dest='BFfilter_keep_complexity', action='store', nargs='?',
                    help='Upper and lower limits for filtering with sequence complexity (default: 0,2).'
                         'Input "off" to disable filtering by sequence complexity.',
                    default='0,2', metavar='STRING')
parser.add_argument('--BF_dimer', dest='BFfilter_keep_dimer', action='store', nargs='?',
                    help='Upper and lower percentile limits for filtering with dimer structure (default: 0,0.85)'
                         'Input "off" to disable filtering by sequence complexity.',
                    default='0,0.85', metavar='STRING')
parser.add_argument('--BF_tm_method', dest='BFfilter_tm_method', action='store', nargs='?',
                    help='Method of calculating melting temperature (tm_NN/tm_GC; default: tm_NN).',
                    choices=['tm_NN', 'tm_GC'], default='tm_NN', metavar='STRING')
parser.add_argument('--BF_tm_table', dest='BFfilter_tm_table', action='store', nargs='?',
                    help='Table of thermodynamic NN values.'
                         'DNA/DNA hybridizations:'
                         'DNA_NN1: values from Breslauer et al. (1986)'
                         'DNA_NN2: values from Sugimoto et al. (1996)'
                         'DNA_NN3: values from Allawi & SantaLucia (1997)'
                         'DNA_NN4: values from SantaLucia & Hicks (2004)'
                         'For RNA/DNA hybridizations:'
                         'R_DNA_NN1: values from Sugimoto et al. (1995) (default).',
                    choices=['DNA_NN1', 'DNA_NN2', 'DNA_NN3', 'DNA_NN4', 'R_DNA_NN1'],
                    default='R_DNA_NN1', metavar='STRING')
parser.add_argument('--BF_hp_ms', dest='BFfilter_match_score', action='store', nargs='?',
                    help='Match score for calculating hairpin score (default: 1).',
                    default='1', metavar='STRING')
parser.add_argument('--BF_hp_mmp', dest='BFfilter_mismatch_penalty', action='store', nargs='?',
                    help='Mismatch penalty for calculating hairpin score (default: -1).',
                    default='-1', metavar='STRING')
parser.add_argument('--BF_hp_gop', dest='BFfilter_gapopen_penalty', action='store', nargs='?',
                    help='Gap open penalty for calculating hairpin score (default: -5).',
                    default='-5', metavar='STRING')
parser.add_argument('--BF_hp_gep', dest='BFfilter_gapextend_penalty', action='store', nargs='?',
                    help='Gap extension penalty for calculating hairpin score (default: -5).',
                    default='-5', metavar='STRING')
parser.add_argument('--BF_dm_kmer', dest='BFfilter_kmer_length', action='store', nargs='?',
                    help='Length of kmer for calculating dimer score (default: 11).',
                    default='11', metavar='STRING')
parser.add_argument('--BF_dm_min_freq', dest='BFfilter_kmer_minfreq', action='store', nargs='?',
                    help='Minimum value of kmer frequency for calculating dimer score (default: 2).',
                    default='2', metavar='STRING')

args = parser.parse_args()

# required input:
SNP_DF = args.snp_df
OUTPUT_BASE = args.output
REFERENCE = args.reference
THREAD = args.thread
LENGTH = args.length
SHIFT = args.shift

# optional:
BG_DATABASE = args.BGfilter_database
BG_MIN_HITS = args.BGfilter_min_hits

AC_DATABASE = args.ACfilter_database
AC_KEEP_HIT = args.ACfilter_keep_hit
AC_MAP_MODE = args.ACfilter_map_mode
AC_TRIM5 = args.ACfilter_trim5
AC_TRIM3 = args.ACfilter_trim3
AC_FILTER_MODE = args.ACfilter_mode
AC_MIN_MAPQ = args.ACfilter_min_mapq
AC_MAX_NM = args.ACfilter_max_nm
AC_MAX_XO = args.ACfilter_max_xo
AC_MAX_XG = args.ACfilter_max_xg
AC_MAX_XSR = args.ACfilter_max_xsr

TX_DATABASE = args.TXfilter_database
TX_KEEP_HIT = args.TXfilter_keep_hit
TX_NAMES_DMP = args.TXfilter_names_dmp
TX_NODES_DMP = args.TXfilter_nodes_dmp
TX_ACC2TAX = args.TXfilter_acc2tax
TX_TAXA = args.TXfilter_taxa
TX_MIN_EDIT = args.TXfilter_min_edit
TX_MAX_EDIT = args.TXfilter_max_edit

BF_TM_METHOD = args.BFfilter_tm_method
BF_TM_TABLE = args.BFfilter_tm_table
BF_HAIRPIN_MATCH_SCORE = args.BFfilter_match_score
BF_HAIRPIN_MISMATCH_PENALTY = args.BFfilter_mismatch_penalty
BF_HAIRPIN_GAPOPEN_PENALTY = args.BFfilter_gapopen_penalty
BF_HAIRPIN_GAPEXTEND_PENALTY = args.BFfilter_gapextend_penalty
BF_DIMER_KMER = args.BFfilter_kmer_length
BF_DIMER_MIN_FREQ = args.BFfilter_kmer_minfreq
BF_KEEP_GC = args.BFfilter_keep_gc
BF_KEEP_TM = args.BFfilter_keep_tm
BF_KEEP_HAIRPIN = args.BFfilter_keep_hairpin
BF_KEEP_COMPLEXITY = args.BFfilter_keep_complexity
BF_KEEP_DIMER = args.BFfilter_keep_dimer

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, 'eModules')
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

log_filename = f'{OUTPUT_BASE}.filter.log'
if os.path.exists(log_filename):
    os.remove(log_filename)
logging.basicConfig(filename=log_filename, format='%(message)s', level=logging.INFO)

TEMP_PROBE = f'{OUTPUT_BASE}.temp_probe.fasta'
PROCESSING_PROBE = f'{OUTPUT_BASE}.processing_probe.fasta'

get_probes_script = 'Get_probe_from_snp.py'
get_probes_params = ['-d', SNP_DF, '-g', REFERENCE, '-l', LENGTH, '-s', SHIFT, '-o', OUTPUT_BASE]
run_python_script(get_probes_script, PATH_SEARCH_LIST, get_probes_params)

check_file_exists(PROCESSING_PROBE, success=f'Successfully converted SNPs to probe sequences.',
                  fail=f'Failed to convert SNPs to probe sequences.')

if BG_DATABASE:
    print('Start background noise filtering...')
    BGfilter_script = 'Snp_filter_bg.py'
    os.renames(PROCESSING_PROBE, TEMP_PROBE)
    BGfilter_params = ['-f', TEMP_PROBE, '-d', BG_DATABASE, '-t', THREAD, '--min_hits', BG_MIN_HITS, '-o',
                       OUTPUT_BASE]
    run_python_script(BGfilter_script, PATH_SEARCH_LIST, BGfilter_params)

    check_file_exists(PROCESSING_PROBE, success=f'Background noise filtering completed.',
                      fail=f'Background noise filtering failed.')
else:
    print('Turn off background noise filter. Skipping...')

if AC_DATABASE:
    print('Start filtering SNPs in inaccessible regions...')
    ACfilter_script = 'Snp_filter_ac_multiprocessing.py'
    os.renames(PROCESSING_PROBE, TEMP_PROBE)
    ACfilter_params = ['-f', TEMP_PROBE, '-i', AC_DATABASE, '-k', AC_KEEP_HIT, '-t', THREAD,
                       '--mm', AC_MAP_MODE, '--trim5', AC_TRIM5, '--trim3', AC_TRIM3, '--fm', AC_FILTER_MODE,
                       '--mapq',
                       AC_MIN_MAPQ, '--nm', AC_MAX_NM, '--xo', AC_MAX_XO, '--xg', AC_MAX_XG, '--xsr', AC_MAX_XSR,
                       '-o',
                       OUTPUT_BASE]
    run_python_script(ACfilter_script, PATH_SEARCH_LIST, ACfilter_params)

    check_file_exists(PROCESSING_PROBE, success=f'Accessibility filtering completed.',
                      fail=f'Accessibility filtering failed.')

else:
    print('Turn off accessibility filter. Skipping...')

if TX_DATABASE:
    print('Start taxonomic filtering...')
    TXfilter_script = 'Snp_filter_tx.py'
    os.renames(PROCESSING_PROBE, TEMP_PROBE)
    TXfilter_params = ['-f', TEMP_PROBE, '--thread', THREAD, '-d', TX_DATABASE, '-k',
                       TX_KEEP_HIT, '--taxa', TX_TAXA, '--names', TX_NAMES_DMP, '--nodes', TX_NODES_DMP, '-a',
                       TX_ACC2TAX, '--minedit', TX_MIN_EDIT, '--maxedit', TX_MAX_EDIT, '-o', OUTPUT_BASE]
    run_python_script(TXfilter_script, PATH_SEARCH_LIST, TXfilter_params)

    check_file_exists(PROCESSING_PROBE, success=f'Taxonomic filtering completed.',
                      fail=f'Taxonomic filtering failed.')

else:
    print('Turn off taxonomic filter. Skipping...')

if args.BFfilter == 'on':
    print('Start filtering with biophysical/biochemical features...')

    os.renames(PROCESSING_PROBE, TEMP_PROBE)

    get_tag_gc_script = 'Get_tag_gc.py'
    get_tag_tm_script = 'Get_tag_tm.py'
    get_tag_dust_script = 'Get_tag_complexity.py'
    get_tag_hairpin_script = 'Get_tag_hairpin.py'
    get_tag_dimer_script = 'Get_tag_dimer.py'

    tmp_dir = f'{OUTPUT_BASE}_eProbe_temp'
    make_temp_dir(tmp_dir)

    tmp_output_base = os.path.join(tmp_dir, os.path.basename(OUTPUT_BASE))

    if BF_KEEP_GC != 'off':
        gc_script_params = ['-f', TEMP_PROBE, '-t', THREAD, '-o',
                            f'{tmp_output_base}.gc']
        run_python_script(get_tag_gc_script, PATH_SEARCH_LIST, gc_script_params)

    if BF_KEEP_TM != 'off':
        tm_script_params = ['-f', TEMP_PROBE, '-t', THREAD, '-m', BF_TM_METHOD, '--table', BF_TM_TABLE,
                            '-o', f'{tmp_output_base}.tm']
        run_python_script(get_tag_tm_script, PATH_SEARCH_LIST, tm_script_params)

    if BF_KEEP_COMPLEXITY != 'off':
        complexity_script_params = ['-f', TEMP_PROBE, '-t', THREAD, '-o',
                                    f'{tmp_output_base}.complexity']
        run_python_script(get_tag_dust_script, PATH_SEARCH_LIST, complexity_script_params)

    if BF_KEEP_HAIRPIN != 'off':
        hairpin_script_params = ['-f', TEMP_PROBE, '-t', THREAD,
                                 '--mp', BF_HAIRPIN_MATCH_SCORE, '--mmp', BF_HAIRPIN_MISMATCH_PENALTY, '--gp',
                                 BF_HAIRPIN_GAPOPEN_PENALTY, '--gpe', BF_HAIRPIN_GAPEXTEND_PENALTY,
                                 '-o', f'{tmp_output_base}.hairpin']
        run_python_script(get_tag_hairpin_script, PATH_SEARCH_LIST, hairpin_script_params)

    if BF_KEEP_DIMER != 'off':
        dimer_script_params = ['-f', TEMP_PROBE, '-k', BF_DIMER_KMER, '--min_freq',
                               BF_DIMER_MIN_FREQ, '-t', THREAD, '-o',
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

    else:
        raise ValueError("No tags were selected. Please at least select one tag to filter with.")

    tags_filter_script = 'Tags_filter.py'
    tags_filter_script_params = ['-m', f'{tmp_output_base}.merged_tags.tsv', '--tm', BF_KEEP_TM, '--gc', BF_KEEP_GC,
                                 '--hairpin', BF_KEEP_HAIRPIN, '--complexity', BF_KEEP_COMPLEXITY, '--dimer',
                                 BF_KEEP_DIMER,
                                 '-o', f'{OUTPUT_BASE}.filtered_probes']
    run_python_script(tags_filter_script, PATH_SEARCH_LIST, tags_filter_script_params)
    shutil.rmtree(tmp_dir)

    check_file_exists(f'{OUTPUT_BASE}.filtered_probes.tsv', success=f'Biophysical/biochemical filtering completed.',
                      fail=f'Biophysical/biochemical filtering failed.')

else:
    print('Turn off tag filter. Skipping...')

if os.path.exists(PROCESSING_PROBE):
    os.renames(PROCESSING_PROBE, TEMP_PROBE)
    filtered_snp_df = probe_to_SNP(TEMP_PROBE)
    filtered_snp_df.to_csv(OUTPUT_BASE + '.filtered_probes.tsv', sep='\t', index=False)

if os.path.exists(TEMP_PROBE):
    os.remove(TEMP_PROBE)
