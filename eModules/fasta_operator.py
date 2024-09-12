from Bio import SeqIO
import collections
import os
import pandas as pd
import shutil
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# read fasta into dict
def read_fasta(fasta_file):
    fasta_dict = collections.OrderedDict()
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    for record in records:
        fasta_dict[record.id] = str(record.seq).replace('\n', '')
    return fasta_dict


# write dict into fasta
def write_fasta(fasta_dict, output_name):
    records = []
    for seq_id, sequence in fasta_dict.items():
        record = SeqRecord(Seq(sequence), id=seq_id)
        records.append(record)

    with open(output_name, 'w') as output_file:
        for record in records:
            output_file.write(f'>{record.id}\n{record.seq.strip()}\n')


# get length of each chromosome
def get_chr_size(fasta_dict, out_base):
    with open(f'{out_base}_chr_len.tsv', 'w') as out:
        for header, seq in fasta_dict.items():
            out.write(header + '\t' + str(len(seq)) + '\n')
    out.close()
    return f'{out_base}_chr_len.tsv'


# read length of each chromosome
def read_chr_size(chr_size_file):
    genome_size_dict = {}
    with open(chr_size_file, 'r') as file:
        for line in file:
            chrom_id, chrom_len = line.strip().split('\t')
            genome_size_dict[chrom_id] = int(chrom_len)
    return genome_size_dict


# split a fasta into multi subset and store in a temp dir
def split_fasta(input_fasta, num_splits, prefix):
    tmp_dir = f'{prefix}_eProbe_temp'
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir)

    records = list(SeqIO.parse(input_fasta, 'fasta'))

    num_records_per_split = len(records) // num_splits
    if len(records) % num_splits != 0:
        num_records_per_split += 1

    output_files = []
    for i in range(num_splits):
        start = i * num_records_per_split
        end = min((i + 1) * num_records_per_split, len(records))
        output_file = os.path.join(tmp_dir, f'{os.path.basename(prefix)}.split{i + 1}.fasta')
        output_files.append(output_file)
        with open(output_file, 'w') as out_handle:
            SeqIO.write(records[start:end], out_handle, 'fasta')

    return output_files


# split a fasta into multi dicts
def split_fasta_into_dicts(fasta_file, num_dicts):
    fasta_dicts = [{} for _ in range(num_dicts)]
    current_dict = 0

    with open(fasta_file, 'r') as file:
        fasta_id = None
        fasta_sequence = ''

        for line in file:
            line = line.strip()

            if line.startswith('>'):
                if fasta_id is not None:
                    fasta_dicts[current_dict][fasta_id] = fasta_sequence

                fasta_id = line[1:]
                fasta_sequence = ''
                current_dict = (current_dict + 1) % num_dicts
            else:
                fasta_sequence += line

        if fasta_id is not None and fasta_sequence:
            fasta_dicts[current_dict][fasta_id] = fasta_sequence

    return fasta_dicts


def probe_header_parser(header):
    chr = header.split(':')[0]
    start, end = header.split(':')[1].split('_')[0].split('-')
    pos, type_, ref, alt = header.split(':')[1].split('_')[1:5]
    return [chr, start, end, pos, type_, ref, alt]


def probe_to_SNP(fasta):
    columns = ['chr', 'start', 'end', 'pos', 'type', 'ref', 'alt']
    probe_data = []
    fasta_dict = read_fasta(fasta)
    for header in fasta_dict.keys():
        probe_data.append(dict(zip(columns, probe_header_parser(header))))
    df = pd.DataFrame(probe_data, columns=columns)
    return df
