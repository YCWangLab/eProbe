from Bio.Align import PairwiseAligner
from Bio import SeqIO
import subprocess
import argparse
import collections
import multiprocessing

usage = """python Get_Tag_Hairpin.py -p pre_probes.fasta -s 48 -t 8 -o output"""

parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-p", "--pre_probes", dest="pre_probes", action="store", nargs='?',
                    help="Pre_probes sequence in fasta", metavar="FILE")
parser.add_argument("-s", "--split", dest="split", action="store", nargs='?', help="Number of splitting files",
                    metavar="STRING")
parser.add_argument("-t", "--threads", dest="threads", action="store", nargs='?', help="Number of threads",
                    metavar="STRING")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output SNPs dataframe with hairpin label", metavar="STRING")
parser.add_argument("--mp", dest="match_p", action="store", nargs='?', default=1, help="Match score for PairwiseAlign",
                    metavar="STRING")
parser.add_argument("--mmp", dest="match_p", action="store", nargs='?', default=-1,
                    help="Mismatch score for PairwiseAlign", metavar="STRING")
parser.add_argument("--gp", dest="gap_p", action="store", nargs='?', default=-5, help="Gap score for PairwiseAlign",
                    metavar="STRING")
parser.add_argument("--egp", dest="extend_gap_p", action="store", nargs='?', default=-5,
                    help="Extend gap score for PairwiseAlign", metavar="STRING")
args = parser.parse_args()

Probes_file = args.pre_probes
Thread = int(args.threads)
Output = args.output
Split_num = int(args.split)
# PairwiseAligner parameters
match_s = args.match_p
mismatch_s = args.mismatch_p
gap_s = args.gap_p
extend_gap_s = args.extend_gap_p


# function: split fasta
def split_fasta(fasta_file, num_files):
    with open(fasta_file, 'r') as f:
        fasta_lines = f.readlines()

    # number of files to split
    num_seqs = len([line for line in fasta_lines if line.startswith('>')])
    seqs_per_file = num_seqs // num_files

    # split
    output_files = []
    curr_file = None
    curr_count = 0
    for line in fasta_lines:
        if line.startswith('>'):
            if curr_count % seqs_per_file == 0:
                # 0 or num_seqs == seqs_per_file open a new file
                if curr_file:
                    curr_file.close()
                file_num = curr_count // seqs_per_file + 1
                filename = f"{fasta_file}.split{file_num}.fasta"
                curr_file = open(filename, 'w')
                output_files.append(filename)
            curr_count += 1
        curr_file.write(line)
    curr_file.close()

    return output_files


# function: read reference fasta
def read_fasta(fasta_file):
    fasta_dict = collections.OrderedDict()
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict


def hairpin_PairwiseAligner_score(probes, match=1, mismatch=-1, gap=-5, extend_gap=-5):
    hairpin_score = collections.OrderedDict()
    aligner = PairwiseAligner(mode='local', score="blastn")
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.gap_score = gap
    aligner.extend_gap_score = extend_gap
    for record in probes:
        alignments = aligner.align(record.seq, record.seq, strand="-")
        if len(alignments):
            best_alignment = max(alignments, key=lambda x: (x.score, len(x.target)))
            alignment_sim = str(best_alignment).split("\n")[1].strip()
            gaps = alignment_sim.count("-")
            identities = alignment_sim.count("|")
            mismatches = alignment_sim.count(".")
            score = 3 * int(identities) - 2 * int(mismatches) - 5 * int(gaps)
            # default filter score <= 30, length <= 20 (strict)
            hairpin_score[record.id] = {"score": score}
    return hairpin_score


def process_probe_file(probe_file):
    probes_record = list(SeqIO.parse(probe_file, "fasta"))
    hairpin = hairpin_PairwiseAligner_score(probes_record, match_s, mismatch_s, gap_s, extend_gap_s)
    SNPs_list = []
    for header, info in hairpin.items():
        chr = header.split("_")[0].split(":")[0]
        pos = header.split("_")[1]
        type = header.split("_")[2]
        Hairpin_score = str(info["score"])
        SNPs_list.append([chr, pos, type, Hairpin_score])


if __name__ == '__main__':
    with multiprocessing.Pool(processes=Thread) as pool:
        # 1 Split fasta
        split_files = split_fasta(Probes_file, Split_num)
        # 2 Making labeled SNPs list
        All_SNPs_list = []
        results = [pool.apply_async(process_probe_file, args=(probe_file,)) for probe_file in split_files]
        for result in results:
            All_SNPs_list = All_SNPs_list + result.get()
        SNPs_df = sorted(All_SNPs_list, key=lambda x: (x[0], int(x[1])))
        for file in split_files:
            rm_temp_cmd = "rm %s" % file
            subprocess.call(rm_temp_cmd, stdout=subprocess.PIPE, shell=True)
        with open(Output, "w") as out:
            out.write("\t".join(["chr", "pos", "type", "hairpin"]) + "\n")
            for SNP in SNPs_df:
                out.write("\t".join(SNP) + "\n")
            out.close()
