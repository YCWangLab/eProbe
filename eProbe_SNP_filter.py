import os
import subprocess
import argparse
from eModules.pipe_controler import *
from eModules.fasta_operator import *
import logging
import pandas as pd
import shutil
from Bio import SeqIO
import collections

usage = """python eProbe_SNP_filter.py -f snp.tsv -r ref.fasta -t threads -o output"""
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
parser.add_argument("--probe_length", dest="probe_length", action="store", nargs='?',
                    help="Length of probes (odd number recommended).",
                    default="81", metavar="STRING")
parser.add_argument("--probe_shift", dest="probe_shift", action="store", nargs='?',
                    help="Position shift (default: 0) fo method pos. - for shifting left "
                         "and + for shifting right.",
                    default="0", metavar="STRING")
parser.add_argument("--BGfilter", dest="BGfilter", action="store", nargs='?',
                    help="Enable background filter to filter out probes overlapping "
                         "with unwanted background noise genomes.",
                    choices=["on", "off"], default="on", metavar="STRING")
parser.add_argument("--BG_db", dest="BGfilter_database", action="store", nargs='?',
                    help="Kraken2 databases (separated with comma) for background filtering.",
                    metavar="STRING")
parser.add_argument("--BG_min_hits", dest="BGfilter_min_hits", action="store", nargs='?',
                    help="Minimum number of non-overlapping kmer hits (separate kmer stacks)"
                         "to filter out (default: 2).",
                    default="2", metavar="STRING")
parser.add_argument("--ACfilter", dest="ACfilter", action="store", nargs='?',
                    help="Enable accessibility filter to filter out probes that can not be"
                         "reliably aligned to given genomes",
                    choices=["on", "off"], default="on", metavar="STRING")
parser.add_argument("--AC_db", dest="ACfilter_database", action="store", nargs='?',
                    help="Bowtie2 databases (separated with comma) for accessibility filtering.",
                    metavar="FILE")
parser.add_argument("--AC_keep_hit", dest="ACfilter_keep_hit", action="store", nargs='?',
                    help="The number of hits to be kept in bowtie2 mapping in accessibility "
                         "filtering (default: 100).",
                    default='100', metavar="STRING")
parser.add_argument("--AC_map_mode", dest="ACfilter_map_mode", action="store", nargs='?',
                    help="Bowtie2 local mapping mode in accessibility filtering. "
                         "(very-fast/fast/sensitive/very-sensitive, default: sensitive)",
                    default="sensitive", metavar="STRING")
parser.add_argument("--AC_trim5", dest="ACfilter_trim5", action="store", nargs='?',
                    help="Trim bases from 5'/left end of reads (default: 0).",
                    default="0", metavar="STRING")
parser.add_argument("--AC_trim3", dest="ACfilter_trim3", action="store", nargs='?',
                    help="Trim bases from 3'/right end of reads (default: 0).",
                    default="0", metavar="STRING")
parser.add_argument("--AC_filter_mode", dest="ACfilter_mode", action="store", nargs='?',
                    help="Accessibility filtering mode (strict/moderate/custom, default: strict)."
                         "If custom is enabled, custom filter mode parameters will be used.",
                    default="strict", metavar="STRING")
parser.add_argument("--AC_min_mapq", dest="ACfilter_min_mapq", action="store", nargs='?',
                    help="[custom filter mode parameter] Minimum mapping quality (default: 30).",
                    default="30", metavar="STRING")
parser.add_argument("--AC_max_nm", dest="ACfilter_max_nm", action="store", nargs='?',
                    help="[custom filter mode parameter] Maximum edit distances of alignment (default: 3). ",
                    default="3", metavar="STRING")
parser.add_argument("--AC_max_xo", dest="ACfilter_max_xo", action="store", nargs='?',
                    help="[custom filter mode parameter] Maximum gap number of alignment (default: 0).",
                    default="0", metavar="STRING")
parser.add_argument("--AC_max_xg", dest="ACfilter_max_xg", action="store", nargs='?',
                    help="[custom filter mode parameter] Maximum total gap length of alignment (default: 0).",
                    default="0", metavar="STRING")
parser.add_argument("--AC_max_xsr", dest="ACfilter_max_xsr", action="store", nargs='?',
                    help="[custom filter mode parameter] Maximum secondary alignment score ratio "
                         "compared to best alignment (default: 0.8).",
                    default="0.8", metavar="STRING")
parser.add_argument("--TXfilter", dest="TXfilter", action="store", nargs='?',
                    help="Enable taxonomic filter to filter out probes that can not be successfully"
                         "assigned to target species (default: on).",
                    choices=["on", "off"], default="on", metavar="STRING")
parser.add_argument("--TX_db", dest="TXfilter_database", action="store", nargs='?',
                    help="Bowtie2 databases (separated with comma) for taxonomic filtering.",
                    metavar="STRING")
parser.add_argument("--TX_keep_hit", dest="TXfilter_keep_hit", action="store", nargs='?',
                    help="The number of hits to be kept in bowtie2 mapping (default: 100).",
                    default='100', metavar="STRING")
parser.add_argument("--TX_taxa", dest="TXfilter_taxa", action="store", nargs='?',
                    help="Taxon number of target species/taxa to be extracted from taxonomic identification",
                    metavar="STRING")
parser.add_argument("--TX_names_dmp", dest="TXfilter_names_dmp", action="store", nargs='?',
                    help="NCBI taxonomy names.dmp", metavar="FILE")
parser.add_argument("--TX_nodes_dmp", dest="TXfilter_nodes_dmp", action="store", nargs='?',
                    help="NCBI taxonomy nodes.dmp", metavar="FILE")
parser.add_argument("--TX_acc2tax", dest="TXfilter_acc2tax", action="store", nargs='?',
                    help="Accession2taxid.txt, the match file of genome accession and taxid.", metavar="FILE")
parser.add_argument("--TX_minedit", dest="TXfilter_min_edit", action="store", nargs='?',
                    help="Minimum edit distance to assign a sequence in ngsLCA (default: 0).",
                    default="0", metavar="STRING")
parser.add_argument("--TX_maxedit", dest="TXfilter_max_edit", action="store", nargs='?',
                    help="Maximum edit distance to assign a sequence in ngsLCA. All sequences have edit distances "
                         "with reference genome lower than this value are regarded as the same distance (default: 2).",
                    default="2", metavar="STRING")
parser.add_argument("--BFfilter", dest="BFfilter", action="store", nargs='?',
                    help="Enable biophysical filter to filter out probes with undesirable biophysical properties",
                    choices=["on", "off"], default="on", metavar="STRING")
parser.add_argument("--BF_Tm_method", dest="BFfilter_Tm_method", action="store", nargs='?',
                    help="Method of calculating melting temperature (Tm_NN/Tm_GC; default: Tm_NN).",
                    choices=["Tm_NN", "Tm_GC"], default="Tm_NN", metavar="STRING")
parser.add_argument("--BF_hp_ms", dest="BFfilter_match_score", action="store", nargs='?',
                    help="Match score for calculating hairpin score (default: 1).",
                    default="1", metavar="STRING")
parser.add_argument("--BF_hp_mmp", dest="BFfilter_mismatch_penalty", action="store", nargs='?',
                    help="Mismatch penalty for calculating hairpin score (default: -1).",
                    default="-1", metavar="STRING")
parser.add_argument("--BF_hp_gop", dest="BFfilter_gapopen_penalty", action="store", nargs='?',
                    help="Gap open penalty for calculating hairpin score (default: -5).",
                    default="-5", metavar="STRING")
parser.add_argument("--BF_hp_gep", dest="BFfilter_gapextend_penalty", action="store", nargs='?',
                    help="Gap extension penalty for calculating hairpin score (default: -5).",
                    default="-5", metavar="STRING")
parser.add_argument("--BF_dm_kmer", dest="BFfilter_kmer_length", action="store", nargs='?',
                    help="Length of kmer for calculating dimer score (default: 11).",
                    default="11", metavar="STRING")
parser.add_argument("--BF_dm_min_freq", dest="BFfilter_kmer_minfreq", action="store", nargs='?',
                    help="Minimum value of kmer frequency for calculating dimer score (default: 2).",
                    default="2", metavar="STRING")
parser.add_argument("--BF_Tm", dest="BFfilter_keep_tm", action="store", nargs='?',
                    help="Upper and lower limits for filtering with melting temperature (default: 68,78).",
                    default="60,80", metavar="STRING")
parser.add_argument("--BF_GC", dest="BFfilter_keep_gc", action="store", nargs='?',
                    help="Upper and lower limits for filtering with GC content (default: 40,60).",
                    default="40,60", metavar="STRING")
parser.add_argument("--BF_hairpin", dest="BFfilter_keep_hairpin", action="store", nargs='?',
                    help="Upper and lower limits for filtering with hairpin score (default: 0,30).",
                    default="0,30", metavar="STRING")
parser.add_argument("--BF_dust", dest="BFfilter_keep_dust", action="store", nargs='?',
                    help="Upper and lower limits for filtering with sequence complexity (default: 0.2).",
                    default="0.2", metavar="STRING")
parser.add_argument("--BF_dimer", dest="BFfilter_keep_dimer", action="store", nargs='?',
                    help="Upper and lower percentile limits for filtering with dimer structure (default: 0,0.85)",
                    default="0,0.85", metavar="STRING")

args = parser.parse_args()

# common input:
SNP_DF = args.snp_df
OUTPUT_PREFIX = args.output
REFERENCE = args.reference
THREAD = args.thread
PROBE_LENGTH = args.probe_length
PROBE_SHIFT = args.probe_shift

# optional:
# Background_probe_filter
BG_DATABASE = args.BGfilter_database
BG_MIN_HITS = args.BGfilter_min_hits

# Accessibility_probe_filter
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

# Taxonomic_probe_filter.py
TX_DATABASE = args.TXfilter_database
TX_KEEP_HIT = args.TXfilter_keep_hit
TX_NAMES_DMP = args.TXfilter_names_dmp
TX_NODES_DMP = args.TXfilter_nodes_dmp
TX_ACC2TAX = args.TXfilter_acc2tax
TX_TAXA = args.TXfilter_taxa
TX_MIN_EDIT = args.TXfilter_min_edit
TX_MAX_EDIT = args.TXfilter_max_edit

# 5 Biophysical_filter.py
BF_TM_METHOD = args.BFfilter_Tm_method
BF_HAIRPIN_MATCH_SCORE = args.BFfilter_match_score
BF_HAIRPIN_MISMATCH_PENALTY = args.BFfilter_mismatch_penalty
BF_HAIRPIN_GAPOPEN_PENALTY = args.BFfilter_gapopen_penalty
BF_HAIRPIN_GAPEXTEND_PENALTY = args.BFfilter_gapextend_penalty
BF_DIMER_KMER = args.BFfilter_kmer_length
BF_DIMER_MIN_FREQ = args.BFfilter_kmer_minfreq
BF_KEEP_GC = args.BFfilter_keep_gc
BF_KEEP_TM = args.BFfilter_keep_tm
BF_KEEP_HAIRPIN = args.BFfilter_keep_hairpin
BF_KEEP_DUST = args.BFfilter_keep_dust
BF_KEEP_DIMER = args.BFfilter_keep_dimer

# make script searching list
PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
MODULES_PATH = os.path.join(PATH_OF_THIS_SCRIPT, "eModules")
PATH_SEARCH_LIST = [PATH_OF_THIS_SCRIPT, MODULES_PATH]

# make log file
logging.basicConfig(filename=f"{OUTPUT_PREFIX}.filter.log", level=logging.INFO)

# generate chr_len.tsv
generate_chr_length_tsv(read_fasta(REFERENCE), OUTPUT_PREFIX)

# 1 run Get_probes_from_SNPs.py
get_probes_script_name = "Get_probes_from_SNPs.py"
get_probes_params = ["-d", SNP_DF, "-g", REFERENCE, "-l", PROBE_LENGTH, "-s", PROBE_SHIFT, "-o", OUTPUT_PREFIX]
run_python_script(get_probes_script_name, PATH_SEARCH_LIST, get_probes_params)

# 2 run Background_probe_filter.py
if args.BGfilter == "on":
    BGfilter_script_name = "Background_probe_filter.py"
    os.renames(f"{OUTPUT_PREFIX}.processing_probe.fasta", f"{OUTPUT_PREFIX}.temp_probe.fasta")
    BGfilter_params = ["-f", f"{OUTPUT_PREFIX}.temp_probe.fasta", "-d", BG_DATABASE, "-t", THREAD, "--min_hits",
                       BG_MIN_HITS]
    run_python_script(BGfilter_script_name, PATH_SEARCH_LIST, BGfilter_params)
    os.remove(f"{OUTPUT_PREFIX}.temp_probe.fasta")

# 3 run Accessibility_probe_filter.py
if args.ACfilter == "on":
    ACfilter_script_name = "Accessibility_filter_multiprocessing.py"
    os.renames(f"{OUTPUT_PREFIX}.processing_probe.fasta", f"{OUTPUT_PREFIX}.temp_probe.fasta")
    ACfilter_params = ["-f", f"{OUTPUT_PREFIX}.temp_probe.fasta", "-i", AC_DATABASE, "-k", AC_KEEP_HIT, "-t", THREAD,
                       "--mm", AC_MAP_MODE, "--trim5", AC_TRIM5, "--trim3", AC_TRIM3, "--fm", AC_FILTER_MODE, "--mapq",
                       AC_MIN_MAPQ, "--nm", AC_MAX_NM, "--xo", AC_MAX_XO, "--xg", AC_MAX_XG, "--xsr", AC_MAX_XSR, "-o",
                       OUTPUT_PREFIX]
    run_python_script(ACfilter_script_name, PATH_SEARCH_LIST, ACfilter_params)
    os.remove(f"{OUTPUT_PREFIX}.temp_probe.fasta")

# 4 run Taxonomic_probe_filter.py
if args.TXfilter == "on":
    TXfilter_script_name = "Taxonomic_probe_filter.py"
    os.renames(f"{OUTPUT_PREFIX}.processing_probe.fasta", f"{OUTPUT_PREFIX}.temp_probe.fasta")
    TXfilter_params = ["-f", f"{OUTPUT_PREFIX}.temp_probe.fasta", "--thread", THREAD, "-d", TX_DATABASE, "-k",
                       TX_KEEP_HIT, "--taxa", TX_TAXA, "--names", TX_NAMES_DMP, "--nodes", TX_NODES_DMP, "-a",
                       TX_ACC2TAX, "--minedit", TX_MIN_EDIT, "--maxedit", TX_MAX_EDIT]
    run_python_script(TXfilter_script_name, PATH_SEARCH_LIST, TXfilter_params)
    os.remove(f"{OUTPUT_PREFIX}.temp_probe.fasta")

# 5 run Biophysical_filter.py
if args.BFfilter == "on":
    get_tag_gc_script = "Get_tag_GC.py"
    get_tag_tm_script = "Get_tag_Tm.py"
    get_tag_dust_script = "Get_tag_DUST.py"
    get_tag_hairpin_script = "Get_tag_hairpin.py"
    get_tag_dimer_script = "Get_tag_dimer.py"

    gc_script_params = ["-f", f"{OUTPUT_PREFIX}.processing_probe.fasta", "-t", THREAD, "-o", OUTPUT_PREFIX + "_gc.tsv"]
    run_python_script(get_tag_gc_script, PATH_SEARCH_LIST, gc_script_params)

    tm_script_params = ["-f", f"{OUTPUT_PREFIX}.processing_probe.fasta", "-t", THREAD, "-m", BF_TM_METHOD, "-o",
                        OUTPUT_PREFIX + "_tm.tsv"]
    run_python_script(get_tag_tm_script, PATH_SEARCH_LIST, tm_script_params)

    dust_script_params = ["-f", f"{OUTPUT_PREFIX}.processing_probe.fasta", "-t", THREAD, "-o",
                          OUTPUT_PREFIX + "_dust.tsv"]
    run_python_script(get_tag_dust_script, PATH_SEARCH_LIST, dust_script_params)

    hairpin_script_params = ["-f", f"{OUTPUT_PREFIX}.processing_probe.fasta", "-t", THREAD,
                             "--mp", BF_HAIRPIN_MATCH_SCORE, "--mmp", BF_HAIRPIN_MISMATCH_PENALTY, "--gp",
                             BF_HAIRPIN_GAPOPEN_PENALTY, "--gpe", BF_HAIRPIN_GAPEXTEND_PENALTY,
                             "-o", OUTPUT_PREFIX + "_hairpin.tsv"]
    run_python_script(get_tag_hairpin_script, PATH_SEARCH_LIST, hairpin_script_params)

    dimer_script_params = ["-f", f"{OUTPUT_PREFIX}.processing_probe.fasta", "-k", BF_DIMER_KMER, "--min_freq",
                           BF_DIMER_MIN_FREQ, "-t", THREAD, "-o", OUTPUT_PREFIX + "_dimer.tsv"]
    run_python_script(get_tag_dimer_script, PATH_SEARCH_LIST, dimer_script_params)

    # merge
    merge_tags_script = "Merge_Tags.py"
    merge_tags_script_params = ["-f", f"{OUTPUT_PREFIX}_gc.tsv", f"{OUTPUT_PREFIX}_tm.tsv", f"{OUTPUT_PREFIX}_dust.tsv",
                                f"{OUTPUT_PREFIX}_hairpin.tsv", f"{OUTPUT_PREFIX}_dimer.tsv", "-o", OUTPUT_PREFIX]
    run_python_script(merge_tags_script, PATH_SEARCH_LIST, merge_tags_script_params)

    # filter
    tags_filter_script = "Tags_filter.py"
    tags_filter_script_params = ["-m", f"{OUTPUT_PREFIX}_merged_tags.tsv", "--Tm", BF_KEEP_TM, "--GC", BF_KEEP_GC,
                                 "--hairpin", BF_KEEP_HAIRPIN, "--dust", BF_KEEP_DUST, "--dimer", BF_KEEP_DIMER, "-o",
                                 f"{OUTPUT_PREFIX}_filtered_SNP"]
    run_python_script(tags_filter_script, PATH_SEARCH_LIST, tags_filter_script_params)

    # remove temp tag files
    tags_files = [f"{OUTPUT_PREFIX}_gc.tsv", f"{OUTPUT_PREFIX}_tm.tsv", f"{OUTPUT_PREFIX}_dust.tsv",
                  f"{OUTPUT_PREFIX}_hairpin.tsv", f"{OUTPUT_PREFIX}_dimer.tsv", f"{OUTPUT_PREFIX}_merged_tags.tsv"]
    delete_files(tags_files)
else:
    filtered_snp_df = probes_to_SNPs(f"{OUTPUT_PREFIX}.processing_probe.fasta")
    filtered_snp_df.to_csv(OUTPUT_PREFIX + "filtered_SNP.tsv", sep='\t', index=False)

os.remove(f"{OUTPUT_PREFIX}.processing_probe.fasta")
