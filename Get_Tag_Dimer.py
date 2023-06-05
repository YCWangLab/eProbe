from Bio import SeqIO
import collections
import subprocess
import argparse
import multiprocessing

usage = """python Get_Tag_Dimer.py -p pre_probes.fasta -s split_num -m method -t threads -o output.txt"""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument("-p", "--pre_probes", dest="pre_probes", action="store", nargs='?',
                    help="Pre_probes sequence in fasta", metavar="FILE")
parser.add_argument("-s", "--split", dest="split", action="store", nargs='?', help="Number of splitting files",
                    metavar="STRING")
parser.add_argument("-t", "--threads", dest="threads", action="store", nargs='?', help="Number of threads",
                    metavar="STRING")
parser.add_argument("-k", "--kmer", dest="kmer_length", action="store", default=9, nargs='?', help="Length of kmer",
                    metavar="STRING")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output SNPs dataframe with dimer label", metavar="FILE")

args = parser.parse_args()

Probes_file = args.pre_probes
Split_num = int(args.split)
Thread = int(args.threads)
Output = args.output
kmer_len = int(args.kmer_length)


# Function: split fasta file
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


# Kmer method functions


# Function: reverse complement
def reverse_complement(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    RCseq_dict = collections.OrderedDict()
    for record in sequences:
        id = record.id
        seq = record.seq.reverse_complement()
        RCseq_dict[id] = str(seq)
    return RCseq_dict


# Function: generate kmer for a sequence
def process_sequence(sequence, k, kmer_counts):
    n = len(sequence)
    for i in range(n - k + 1):
        kmer = sequence[i:i + k]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1


# Function: kmer dict generator
def generate_kmers(fasta_dict, k):
    kmer_counts = {}
    for seq in fasta_dict.values():
        process_sequence(seq, k, kmer_counts)
    return kmer_counts


# Function: calculate dimer score for a sequence
def search_kmers_in_dict(sequence, k, dictionary):
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    total_sum = 0

    for kmer in kmers:
        if kmer in dictionary:
            total_sum += dictionary[kmer]
    return total_sum


# Function: calculate dimer score for a fasta file
def dimer_kmer_score(probes_file, k, kmer_dict):
    kmer_score = {}
    sequences = list(SeqIO.parse(probes_file, "fasta"))
    for record in sequences:
        id = record.id
        score = search_kmers_in_dict(str(record.seq), k, kmer_dict)
        kmer_score[id] = score
    return kmer_score


# Function: merge dict
def merge_dicts(*dicts):
    merged_dict = {}
    for dictionary in dicts:
        for key, value in dictionary.items():
            if key in merged_dict:
                merged_dict[key] += value
            else:
                merged_dict[key] = value
    return merged_dict


def process_probe_file(probe_file, kmer_dict=None):
    dimer_dict = dimer_kmer_score(probe_file, kmer_len, kmer_dict)
    SNPs_list = []
    for header, score in dimer_dict.items():
        chr = header.split("_")[0].split(":")[0]
        pos = header.split("_")[1]
        type = header.split("_")[2]
        Dimer_score = str(score)
        SNPs_list.append([chr, pos, type, Dimer_score])
    return SNPs_list


if __name__ == '__main__':
    with multiprocessing.Pool(processes=Thread) as pool:
        # 1 Split fasta
        split_files = split_fasta(Probes_file, Split_num)
        # 2 Making kmer dict
        results = []
        for file in split_files:
            # reverse complement
            rcseq_dict = reverse_complement(file)
            results.append(pool.apply_async(generate_kmers, args=(rcseq_dict, kmer_len)))
        All_kmer_dict = {}
        for result in results:
            All_kmer_dict = merge_dicts(All_kmer_dict, result.get())
        # 3 Making SNPs list
        All_SNPs_list = []
        results = []
        for file in split_files:
            # Making labeled SNPs list
            results.append(pool.apply_async(process_probe_file, args=(file, All_kmer_dict)))
        for result in results:
            All_SNPs_list = All_SNPs_list + result.get()
        SNPs_df = sorted(All_SNPs_list, key=lambda x: (x[0], int(x[1])))
        for file in split_files:
            rm_temp_cmd = "rm %s" % file
            subprocess.call(rm_temp_cmd, stdout=subprocess.PIPE, shell=True)
        # output
        with open(Output, "w") as out:
            out.write("\t".join(["chr", "pos", "type", "dimer"]) + "\n")
            for SNP in SNPs_df:
                out.write("\t".join(SNP) + "\n")
            out.close()
