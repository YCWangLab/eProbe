from Bio.SeqUtils import MeltingTemp as Tm
from Bio.SeqUtils import GC
import collections
from Bio import SeqIO
import argparse
import multiprocessing

usage = """python Get_Tag_GC_Tm.py -p pre_probes.fasta -m "Tm_NN" -t threads -o output.txt"""

parser = argparse.ArgumentParser(description=usage)

parser.add_argument("-p", "--pre_probes", dest="pre_probes", action="store", nargs='?',
                    help="Pre_probes sequence in fasta", metavar="FILE")
parser.add_argument("-m", "--Tm_method", dest="Tm_method", action="store", nargs='?',
                    help="Method for melting temperature calculation (Tm_GC/Tm_NN)", metavar="STRING")
parser.add_argument("-t", "--threads", dest="threads", action="store", nargs='?', help="Number of threads",
                    metavar="STRING")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output SNPs dataframe with GC and Tm labels", metavar="FILE")

args = parser.parse_args()

Probes_file = args.pre_probes
Tm_method = args.Tm_method
Thread = int(args.threads)
Output = args.output


# function: read reference fasta
def read_fasta(fasta_file):
    fasta_dict = collections.OrderedDict()
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict


# Tm method
if Tm_method == "Tm_NN":
    def Tm_cal(seq):
        # Allawi & SantaLucia (1997)
        t = Tm.Tm_NN(seq, nn_table=Tm.DNA_NN3)
        return t
elif Tm_method == "Tm_GC":
    def Tm_cal(seq):
        t = Tm.Tm_GC(seq)
        return t


def process_probe(header, seq):
    chr = header.split("_")[0].split(":")[0]
    pos = header.split("_")[1]
    type = header.split("_")[2]
    gc = str(round(GC(seq), 4))
    tm = str(round(Tm_cal(seq), 4))
    return [chr, pos, type, gc, tm]


# read probes file
probes = read_fasta(Probes_file)

# Making labeled SNPs data frame
SNPs_list = []

if __name__ == '__main__':
    with multiprocessing.Pool(processes=Thread) as pool:
        results = []
        for header, seq in probes.items():
            results.append(pool.apply_async(process_probe, (header, seq)))
        for result in results:
            SNPs_list.append(result.get())
        SNPs_df = sorted(SNPs_list, key=lambda x: (x[0], int(x[1])))

        with open(Output, "w") as out:
            out.write("\t".join(["chr", "pos", "type", "GC", "Tm"]) + "\n")
            for SNP in SNPs_df:
                out.write("\t".join(SNP) + "\n")
            out.close()
