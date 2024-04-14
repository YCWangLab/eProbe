import subprocess
import pandas as pd
import os
import shutil
import datetime
from pipe_controler import diff_file_line_count


def run_mapping(thread, databse, fasta, output, k):
    bname = os.path.splitext(os.path.basename(databse))[0]
    bowtie2_cmd = f"bowtie2 -f -k {k} --threads {thread} -x {databse} -U {fasta} --no-unal -S {output}.{bname}.sam"
    subprocess.run(bowtie2_cmd, shell=True, check=True)


def filter_header(header_full, header_new, header_subset):
    df1 = pd.read_csv(header_new, header=None)
    df1 = df1[0].unique()
    df1 = [f"SN:{value}" for value in df1]

    header = pd.read_csv(header_full, header=None, sep='\t')
    header = header[header[1].isin(df1)]

    header.to_csv(header_subset, header=False, index=False, sep='\t')


def process_sam_files(bname, bDB, thread):
    sam_file = f"{bname}.{bDB}.sam"
    header_subset_1_file = f"{bname}.{bDB}.header_subset.1.txt"
    header_subset_2_file = f"{bname}.{bDB}.header_subset.2.txt"
    header_subset_3_file = f"{bname}.{bDB}.header_subset.3.txt"
    alignment_file = f"{bname}.{bDB}.alignment.txt"
    header_new_file = f"{bname}.{bDB}.header_new.txt"
    bam_file = f"{bname}.{bDB}.bam"

    subprocess.run(["grep", "@HD", sam_file], stdout=open(header_subset_1_file, "w"))
    subprocess.run(["grep", "@SQ", sam_file], stdout=open(header_subset_2_file, "w"))
    subprocess.run(["grep", "@PG", sam_file], stdout=open(header_subset_3_file, "w"))
    subprocess.run(["grep", "-v", "^@", sam_file], stdout=open(alignment_file, "w"))

    awk_command = f"awk '!seen[$0]++' {header_subset_2_file} > {header_subset_2_file}.tmp"
    subprocess.run(awk_command, shell=True, check=True)
    os.remove(header_subset_2_file)
    os.rename(header_subset_2_file + ".tmp", header_subset_2_file)

    extract_cmd = f"cat {alignment_file} | cut -f3 | sort -u | uniq > {header_new_file}"
    subprocess.run(extract_cmd, shell=True, check=True)

    os.rename(header_subset_2_file, header_subset_2_file + ".tmp")
    filter_header(header_subset_2_file + ".tmp", header_new_file, header_subset_2_file)

    subprocess.run(
        f"cat {header_subset_1_file} {header_subset_2_file} {header_subset_3_file} {alignment_file} | samtools view -@ {thread} -b -o {bam_file}",
        shell=True, check=True)

    os.remove(header_subset_1_file)
    os.remove(header_subset_2_file)
    os.remove(header_subset_3_file)
    os.remove(sam_file)


def merge_sort(bam_path, bname, thread):
    bam_files = [f for f in os.listdir(bam_path) if f.endswith('.bam')]
    if len(bam_files) == 1:
        merged_bam = os.path.join(bam_path, bam_files[0])
    else:
        merge_cmd = f"samtools merge -n -@ {thread} {bam_path}/{bname}.merged {' '.join([f'{bam_path}/{f}' for f in bam_files])}"
        subprocess.run(merge_cmd, shell=True, check=True)
        merged_bam = f"{bam_path}/{bname}.merged"

    sort_cmd = f"samtools sort -n -@ {thread} -O bam -o {bname}.sorted.bam {merged_bam}"
    subprocess.run(sort_cmd, shell=True, check=True)

    return f"{bname}.sorted.bam"


def run_ngsLCA(bam, names_dmp, nodes_dmp, acc2tax, minedit, maxedit, output):
    ngsLCA_cmd = f"ngsLCA -names {names_dmp} -nodes {nodes_dmp} -acc2tax {acc2tax} -editdistmin {minedit} -editdistmax {maxedit} -bam {bam} -outnames {output}.min{minedit}max{maxedit}"
    subprocess.run(ngsLCA_cmd, shell=True, check=True)
    return f"{output}.min{minedit}max{maxedit}.lca"


def run_script(script_path, parameters):
    command = ['python', script_path] + parameters
    subprocess.check_output(command, universal_newlines=True)


if __name__ == '__main__':
    import argparse

    usage = """ python Taxonomic_probe_filter.py -f probes.fasta --thread 16 -d ref1,ref2,... --taxa taxid/name --names names.dmp --nodes nodes.dmp -a acc2tax.txt --minedit 0 --maxedit 2 -o output """

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-f", "--fasta", dest="fasta", action="store", nargs='?',
                        help="Input fasta.", metavar="FILE", required=True)
    parser.add_argument("--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Number of thread (default: 1). ", metavar="STRING")
    parser.add_argument("-d", "--database", dest="database", action="store", nargs='?',
                        help="Bowtie2 database(s) (separated with comma).", metavar="STRING", required=True)
    parser.add_argument("-k", dest="keep_hit", action="store", nargs='?', default='100',
                        help="The number of hits to be kept in bowtie2 mapping (default: 100).", metavar="STRING")
    parser.add_argument("--taxa", dest="taxa", action="store", nargs='?',
                        help="Taxon number/latin names of target species/taxa to be extracted from lca",
                        metavar="STRING",
                        required=True)
    parser.add_argument("--names", dest="names", action="store", nargs='?',
                        help="NCBI taxonomy names.dmp", metavar="FILE", required=True)
    parser.add_argument("--nodes", dest="nodes", action="store", nargs='?',
                        help="NCBI taxonomy nodes.dmp", metavar="FILE", required=True)
    parser.add_argument("-a", "--acc2tax", dest="acc2tax", action="store", nargs='?',
                        help="accession2taxid.txt, the match file of genome accession and taxid.", metavar="FILE",
                        required=True)
    parser.add_argument("--minedit", dest="minedit", action="store", nargs='?', default=0,
                        help="Minimum edit distance to assign a sequence in ngsLCA (default: 0).", metavar="STRING",
                        required=True)
    parser.add_argument("--maxedit", dest="maxedit", action="store", nargs='?', default=2,
                        help="Maximum edit distance to assign a sequence in ngsLCA. All sequences have edit distances with reference genome lower than this value are regarded as the same distance (default: 2).",
                        metavar="STRING", required=True)
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)

    args = parser.parse_args()

    tmp_dir = f"{args.output}_eProbe_temp"
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)

    output_path = os.path.join(tmp_dir, args.output)

    for db in args.database.split(","):
        print("Mapping probes to database: %s with bowtie2 ..." % db, datetime.datetime.now())
        run_mapping(int(args.thread), db, args.fasta, output_path, args.keep_hit)
    print("Mapping completed ...", datetime.datetime.now())

    for db in args.database.split(","):
        bname = os.path.splitext(os.path.basename(db))[0]
        print("Preparing input bam file for generated with database: %s ..." % db, datetime.datetime.now())
        process_sam_files(output_path, bname, int(args.thread))
    print("Input bam(s) preparation completed ...", datetime.datetime.now())

    print("Merging and sorting bam ...", datetime.datetime.now())
    prepared_bam = merge_sort(tmp_dir, args.output, int(args.thread))
    print("Merging and sorting completed ...", datetime.datetime.now())

    print("Running taxonomic identification with ngsLCA ...", datetime.datetime.now())
    lca_file = run_ngsLCA(prepared_bam, args.names, args.nodes, args.acc2tax, args.minedit, args.maxedit, args.output)
    print("Taxonomic identification completed ...", datetime.datetime.now())

    print("Extracting probes passing taxonomic identification ...", datetime.datetime.now())
    script_directory = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_directory, 'Extract_seq_from_LCA.py')
    parameter_template = ["-l", lca_file, "-t", args.taxa, "-i", args.fasta, "-f", "fasta",
                          "-o", f"{args.output}.processing_probe.fasta"]
    run_script(script_path, parameter_template)

    before_filtering = str(int(int(subprocess.check_output(["wc", "-l", args.fasta]).split()[0]) / 2))
    after_filtering = str(int(diff_file_line_count(args.fasta, f"{args.output}.processing_probe.fasta") / 2))
    print("Filtered out %s probes from %s probes in taxonomic filtering." % (after_filtering, before_filtering),
          datetime.datetime.now())

    shutil.rmtree(tmp_dir)
    os.remove(prepared_bam)
