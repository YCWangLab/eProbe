import pandas as pd
import vcf


def concat_dfs(file_paths):
    first_df = pd.read_csv(file_paths[0], sep='\t')
    if 'chr' not in first_df.columns or 'pos' not in first_df.columns:
        raise ValueError("The first file must contain 'chr' and 'pos' columns.")

    merged_df = first_df

    for path in file_paths[1:]:
        df = pd.read_csv(path, sep='\t')
        if 'chr' not in df.columns or 'pos' not in df.columns:
            raise ValueError(f"{path} is missing 'chr' or 'pos' columns.")
        merged_df = pd.concat([merged_df, df], ignore_index=True)

    return merged_df


def SNPdf_to_vcf(snp_df_, vcf_reader_, out_vcf):
    vcf_writer = vcf.Writer(open(out_vcf + ".vcf", "w"), template=vcf_reader_)

    for row in snp_df_.itertuples():
        chrom = row.chr
        pos = int(row.pos)

        for record in vcf_reader_.fetch(chrom, pos - 1, pos + 1):
            vcf_writer.write_record(record)

    vcf_writer.close()


if __name__ == '__main__':
    import argparse

    usage = """python Get_VCF_from_SNPs.py -v raw.vcf.gz(*.tbi) -o out.vcf -f df1.txt df2.txt [...]"""
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-f', '--snp_df', dest="snp_df", action="store", nargs='+',
                        help="SNPs dataframe(s).", metavar="FILE", required=True)
    parser.add_argument("-v", "--vcf", dest="vcf", action="store", nargs='?',
                        help="Raw VCF file.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING")
    args = parser.parse_args()

    merged_df = concat_dfs(args.df_files)
    vcf_reader = vcf.Reader(filename=args.vcf, compressed=True)
    SNPdf_to_vcf(merged_df, vcf_reader, args.output)
