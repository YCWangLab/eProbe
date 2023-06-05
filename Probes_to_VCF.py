import pandas as pd
import vcf
import argparse
usage = """python Probes_to_VCF.py -v raw.vcf.gz(*.tbi) -o out.vcf -f df1.txt df2.txt [...]"""
parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-f', '--df_files', dest="df_files", action="store", nargs='+',
                    help="Multiple SNPs dataframes", metavar="FILE")
parser.add_argument("-v", "--vcf", dest="vcf", action="store", nargs='?',
                    help="Raw VCF file", metavar="FILE")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output file name", metavar="FILE")
args = parser.parse_args()
df_files = args.df_files
vcf_file = args.vcf
output = args.output
# function: merge multiple SNPs dataframes
def merge_dataframes(file_paths):
    # reading first data frame
    df = pd.read_csv(file_paths[0], sep='\t', header=0)

    # read and merge the remaining data frame one by one
    for path in file_paths[1:]:
        df2 = pd.read_csv(path, sep='\t', header=0)
        df = pd.merge(df, df2, how='outer')

    return df

# function: output vcf for examination
def SNPdf_to_vcf(snp_df_, vcf_reader_, out_vcf):
    vcf_writer = vcf.Writer(open(out_vcf, "w"), template=vcf_reader_)
    for row in snp_df_.itertuples():
        chrom = row.chr
        pos = int(row.pos)
        for record in vcf_reader_.fetch(chrom, pos - 1, pos):
            vcf_writer.write_record(record)
    vcf_writer.close()

# merged SNPs dataframes
merged_df = merge_dataframes(df_files)
# read vcf
vcf_reader = vcf.Reader(filename=vcf_file, compressed=True)
# output vcf for examination
SNPdf_to_vcf(merged_df, vcf_reader, output)



