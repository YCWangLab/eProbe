import multiprocessing
import pandas as pd


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

                fasta_id = line[1:] # Update the sequence ID.
                fasta_sequence = "" # Reset sequence data.
                current_dict = (current_dict + 1) % num_dicts
            else:
                fasta_sequence += line

        if fasta_id is not None and fasta_sequence:
            fasta_dicts[current_dict][fasta_id] = fasta_sequence

    return fasta_dicts


def cal_dust_score(seq):
    k = 3
    kmers_list = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmers_list.append(kmer)
    kmers_set = list(set(kmers_list))
    Ct_dict = {}
    for i in kmers_set:
        Ct_dict[i] = kmers_list.count(i)
    all_ct = 0
    for ct in Ct_dict.values():
        all_ct = all_ct + int(ct) * (int(ct) - 1) / 2
    return all_ct / (len(seq) - 1)


def process_probe_dict(probe_dict):
    score_list = []
    for header, seq in probe_dict.items():
        chr = header.split(":")[0]
        start, end = header.split(":")[1].split("_")[0].split("-")
        pos, type_, ref, alt = header.split(":")[1].split("_")[1:5]
        tag = str(cal_dust_score(seq))
        score_list.append([chr, start, end, pos, type_, ref, alt, tag])
    return score_list


if __name__ == '__main__':

    import argparse

    usage = """python Get_tag_DUST.py -f probes.fasta -t threads -o output"""
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-f", "--fasta", dest="fasta", action="store", nargs='?',
                        help="Input fasta.", metavar="FILE", required=True)
    parser.add_argument("-t", "--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Number of threads (default: 1).", metavar="STRING")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)
    args = parser.parse_args()

    split_dicts = split_fasta_into_dicts(args.fasta, int(args.thread))

    with multiprocessing.Pool(processes=int(args.thread)) as pool:
        results = [pool.apply_async(process_probe_dict, args=(fa_dict,)) for fa_dict in split_dicts]
        all_SNPs_tag_list = []
        for result in results:
            all_SNPs_tag_list = all_SNPs_tag_list + result.get()
    sorted_SNPs_tag_list = sorted(all_SNPs_tag_list, key=lambda x: (x[0], int(x[1])))
    SNPs_df = pd.DataFrame(sorted_SNPs_tag_list, columns=["chr", "start", "end", "pos", "type", "ref", "alt", "DUST"])
    SNPs_df.to_csv(args.output + ".tsv", sep='\t', index=False)
