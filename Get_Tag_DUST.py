import argparse
import multiprocessing

usage = """python Get_Tag_DUST.py -p pre_probes.fasta -s split_num -t threads -o output.txt"""
parser = argparse.ArgumentParser(description=usage)

parser.add_argument("-p", "--pre_probes", dest="pre_probes", action="store", nargs='?',
                    help="Pre_probes sequence in fasta", metavar="FILE")
parser.add_argument("-s", "--split", dest="split", action="store", nargs='?', help="Number of splitting files",
                    metavar="STRING")
parser.add_argument("-t", "--threads", dest="threads", action="store", nargs='?', help="Number of threads",
                    metavar="STRING")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output SNPs dataframe with complexity label", metavar="FILE")

args = parser.parse_args()

fasta_file = args.pre_probes
Split_num = int(args.split)
Thread = int(args.threads)
Output = args.output


def split_fasta_into_dicts(fasta_file, num_dicts):
    fasta_dicts = [{} for _ in range(num_dicts)]
    current_dict = 0

    with open(fasta_file, 'r') as file:
        fasta_id = None
        fasta_sequence = ""

        for line in file:
            line = line.strip()

            if line.startswith('>'):
                if fasta_id is not None:
                    fasta_dicts[current_dict][fasta_id] = fasta_sequence

                fasta_id = line[1:]
                fasta_sequence = ""
                current_dict = (current_dict + 1) % num_dicts
            else:
                fasta_sequence += line

        if fasta_id is not None and fasta_sequence:
            fasta_dicts[current_dict][fasta_id] = fasta_sequence

    return fasta_dicts


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


# Function calculate dust score
# Formular: A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences (2006)
def cal_dust_score(seq):
    k = 3
    kmers_list = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmers_list.append(kmer)
    kmers_set = list(set(kmers_list))
    # each Ct
    Ct_dict = {}
    for i in kmers_set:
        Ct_dict[i] = kmers_list.count(i)
    # calculate score
    all_ct = 0
    for ct in Ct_dict.values():
        all_ct = all_ct + int(ct) * (int(ct) - 1) / 2
    return all_ct / (len(seq) - 1)


def process_probe_file(probe):
    SNPs_list = []
    for header, seq in probe.items():
        chr = header.split("_")[0].split(":")[0]
        pos = header.split("_")[1]
        type = header.split("_")[2]
        Dust_score = str(cal_dust_score(seq))
        SNPs_list.append([chr, pos, type, Dust_score])
    return SNPs_list


if __name__ == '__main__':
    with multiprocessing.Pool(processes=Thread) as pool:
        # 1 Split fasta
        split_dicts = split_fasta_into_dicts(fasta_file, Thread)
        # 2 Calculate dust score
        results = [pool.apply_async(process_probe_file, args=(fa_dict,)) for fa_dict in split_dicts]
        All_SNPs_list = []
        for result in results:
            All_SNPs_list = All_SNPs_list + result.get()
        SNPs_df = sorted(All_SNPs_list, key=lambda x: (x[0], int(x[1])))
        # 3 Output
        with open(Output, "w") as out:
            out.write("\t".join(["chr", "pos", "type", "DUST"]) + "\n")
            for SNP in SNPs_df:
                out.write("\t".join(SNP) + "\n")
            out.close()
