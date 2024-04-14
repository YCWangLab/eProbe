from Bio import SeqIO
import collections
import os
import pandas as pd

def read_fasta(fasta_file):
    fasta_dict = collections.OrderedDict()
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        fasta_dict[record.id] = str(record.seq.replace("\n", ""))
    return fasta_dict


def generate_chr_length_tsv(fasta_dict, out_pre):
    with open(f"{out_pre}_chr_len.tsv", "w") as out:
        for header, seq in fasta_dict.items():
            out.write(header + "\t" + seq + "\n")
    out.close()


def split_fasta(input_fasta, num_splits, prefix):
    output_dir = f"{prefix}_eProbe_temp"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        # if temp dir exists, rm all files in it
        for file_name in os.listdir(output_dir):
            file_path = os.path.join(output_dir, file_name)
            os.remove(file_path)

    records = list(SeqIO.parse(input_fasta, "fasta"))

    num_records_per_split = len(records) // num_splits
    if len(records) % num_splits != 0:
        num_records_per_split += 1

    output_files = []
    for i in range(num_splits):
        start = i * num_records_per_split
        end = min((i + 1) * num_records_per_split, len(records))
        output_file = os.path.join(output_dir, f'{prefix}.split{i + 1}.fasta')
        output_files.append(output_file)
        with open(output_file, 'w') as out_handle:
            SeqIO.write(records[start:end], out_handle, 'fasta')

    return output_files


def probe_header_parser(header):
    chr = header.split(":")[0]
    start, end = header.split(":")[1].split("_")[0].split("-")
    pos, type_, ref, alt = header.split(":")[1].split("_")[1:5]
    return [chr, start, end, pos, type_, ref, alt]

def probes_to_SNPs(fasta):
    columns = ['chr', 'start', 'end', 'pos', 'type', 'ref', 'alt']
    df = pd.DataFrame(columns=columns)
    fasta_dict = read_fasta(fasta)
    for header in fasta_dict.keys():
        df = df.append(dict(zip(columns, probe_header_parser(header))), ignore_index=True)
    return df
