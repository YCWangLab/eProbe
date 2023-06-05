import collections
import argparse

usage = """python Get_Tag_Taxid.py -k kraken2_output -o output.txt"""

parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-k", "--kraken_output", dest="kraken_output", action="store", nargs='?',
                    help="Kraken2 output file (storing assignment for each sequence)", metavar="FILE")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output SNPs dataframe with taxid label", metavar="STRING")
args = parser.parse_args()

Probes_file = args.kraken_output
Output = args.output


# read header
def taxid_reader(taxid):
    SNPs_dict = collections.OrderedDict()
    with open(taxid, "r") as f:
        for line in f.readlines():
            header = line.split()[1].strip()
            taxid = line.split()[2].strip()
            if header not in SNPs_dict.keys():
                SNPs_dict[header] = []
                SNPs_dict[header].append(taxid)
            else:
                SNPs_dict[header].append(taxid)
    return SNPs_dict


# make SNPs list
def dict2list(SNPs_dict):
    SNPs_list = []
    for header, taxid in SNPs_dict.items():
        chr = header.split("_")[0].split(":")[0]
        pos = header.split("_")[1]
        type = header.split("_")[2].split()[0]
        taxid = ",".join(taxid)
        SNPs_list.append([chr, pos, type, taxid])
    return SNPs_list


# output
SNPs_dict = taxid_reader(Probes_file)
SNPs_list = dict2list(SNPs_dict)
with open(Output, "w") as out:
    out.write("\t".join(["chr", "pos", "type", "taxid"]) + "\n")
    for SNP in SNPs_list:
        out.write("\t".join(SNP) + "\n")
out.close()
