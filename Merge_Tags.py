import collections
import argparse
import pandas as pd

usage = """python Merge_Tags.py -p pre_probes.fasta -g GCTm.txt -a hairpin.txt -d dimer.txt -u dust.txt -t taxid.txt -o output.txt"""
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-p", "--probes", dest="probe_fasta", action="store", nargs='?',
                    help="Probe fasta", metavar="FILE")
parser.add_argument("-g", "--GCTm", dest="label_GCTm", action="store", nargs='?',
                    help="Label file of GC content and melting temperature", metavar="FILE")
parser.add_argument("-a", "--hairpin", dest="label_hairpin", action="store", nargs='?',
                    help="Label file of hairpin structure", metavar="FILE")
parser.add_argument("-d", "--dimer", dest="label_dimer", action="store", nargs='?',
                    help="Label file of dimer structure", metavar="FILE")
parser.add_argument("-u", "--dust", dest="label_dust", action="store", nargs='?',
                    help="Label file of complexity", metavar="FILE")
parser.add_argument("-t", "--taxid", dest="label_taxid", action="store", nargs='?',
                    help="Label file of taxid", metavar="FILE")
parser.add_argument("-o", "--out", dest="output", action="store", nargs='?',
                    help="Merged label file", metavar="FILE")
args = parser.parse_args()

probe_file = args.probe_fasta
label_GCTm = args.label_GCTm
label_hairpin = args.label_hairpin
label_dimer = args.label_dimer
label_dust = args.label_dust
label_taxid = args.label_taxid
output = args.output


def make_SNPs_dict(fasta):
    SNPs_dict = collections.OrderedDict()
    with open(fasta, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                chr = line[1:].split(":")[0].strip()
                pos = line[1:].split("_")[1].strip()
                type = line[1:].split("_")[2].strip()
                SNPs_dict["_".join([chr, pos])] = {"chr": chr, "pos": pos, "type": type, "GC": "N/A", "Tm": "N/A",
                                                   "hairpin": 0,
                                                   "dimer": 0, "DUST": 0, "taxid": 0}

    return SNPs_dict


def read_label(label_file, label, SNPs_dict):
    df = pd.read_csv(label_file, sep='\t', header=0)
    for row in df.itertuples():
        SNPs_dict["_".join([row.chr, str(row.pos)])][label] = str(getattr(row, label))
    return SNPs_dict


# make SNPs_dict
SNPs_dict = make_SNPs_dict(probe_file)
# GC and Tm
SNPs_dict = read_label(label_GCTm, "Tm", SNPs_dict)
SNPs_dict = read_label(label_GCTm, "GC", SNPs_dict)
# hairpin
SNPs_dict = read_label(label_hairpin, "hairpin", SNPs_dict)
# dimer
SNPs_dict = read_label(label_dimer, "dimer", SNPs_dict)
# dust
SNPs_dict = read_label(label_dust, "DUST", SNPs_dict)
# taxid
SNPs_dict = read_label(label_taxid, "taxid", SNPs_dict)
# output
with open(output, "w") as out:
    out.write("\t".join(["chr", "pos", "type", "Tm", "GC", "hairpin", "dimer", "DUST", "taxid"]) + "\n")
    for key, value in SNPs_dict.items():
        out.write("\t".join(map(str,
                                [value["chr"], value["pos"], value["type"], value["Tm"], value["GC"], value["hairpin"],
                                 value["dimer"], value["DUST"], value["taxid"]])) + "\n")
out.close()
