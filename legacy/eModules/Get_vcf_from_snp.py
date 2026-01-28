import pandas as pd
import pysam


def concat_dfs(file_paths):
    first_df = pd.read_csv(file_paths[0], sep='\t')
    if 'chr' not in first_df.columns or 'pos' not in first_df.columns:
        raise ValueError('The first file must contain "chr" and "pos" columns.')

    merged_df = first_df

    for path in file_paths[1:]:
        df = pd.read_csv(path, sep='\t')
        if 'chr' not in df.columns or 'pos' not in df.columns:
            raise ValueError(f'{path} is missing "chr" or "pos" columns.')
        merged_df = pd.concat([merged_df, df], ignore_index=True)

    return merged_df


def SNPdf_to_vcf(snp_df_, vcf_path, out_vcf):
    vcf_in = pysam.VariantFile(vcf_path)
    vcf_out = pysam.VariantFile(out_vcf + '.vcf', 'w', header=vcf_in.header)

    for row in snp_df_.itertuples():
        chrom = row.chr
        pos = int(row.pos)

        for record in vcf_in.fetch(chrom, pos - 1, pos + 1):
            vcf_out.write(record)

    vcf_out.close()
    vcf_in.close()


if __name__ == '__main__':
    import argparse

    usage = '''python Get_vcf_from_snp.py -v raw.vcf.gz(*.tbi) -o out.vcf -f df1.txt df2.txt [...]'''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-f', '--snp_dfs', dest='snp_dfs', action='store', nargs='+',
                        help='Input SNPs data frame(s) (tsv).',
                        metavar='FILE', required=True)
    parser.add_argument('-v', '--vcf', dest='vcf', action='store', nargs='?',
                        help='Gzipped and indexed VCF file (vcf.gz).',
                        metavar='FILE')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='STRING')
    args = parser.parse_args()

    merged_df = concat_dfs(args.snp_dfs)
    SNPdf_to_vcf(merged_df, args.vcf, args.output)
