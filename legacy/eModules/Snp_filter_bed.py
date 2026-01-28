import datetime
import subprocess
from pipe_controler import diff_file_line_count, check_bed


def SNP2BED(snp_df, output):
    awk_command = f'''awk 'BEGIN {{OFS="\\t"}} {{print $1, $2, $2+1, $3, $4, $5}}' {snp_df} > {output}_eProbe_temp.snp'''
    subprocess.run(awk_command, shell=True, check=True)
    return f'{output}_eProbe_temp.snp'


def run_bedtools_slop(bed_file, genome_size, distance, output):
    slop_command = f'''bedtools slop -i {bed_file} -g {genome_size} -l {distance} -r {distance} | awk '($3 - $2) >= 1' > {output}_eProbe_temp.bed'''
    subprocess.run(slop_command, shell=True, check=True)
    return f'{output}_eProbe_temp.bed'


def run_bedtools_keep(snp_bed, bed_file, output):
    keep_command = f'bedtools intersect -a {snp_bed} -b {bed_file} -wa > {output}_eProbe_temp_kept.snp'
    subprocess.run(keep_command, shell=True, check=True)
    return f'{output}_eProbe_temp_kept.snp'


def run_bedtools_remove(snp_bed, bed_file, output):
    rm_command = f'bedtools intersect -v -a {snp_bed} -b {bed_file} -wa > {output}_eProbe_temp_removed.snp'
    subprocess.run(rm_command, shell=True, check=True)
    return f'{output}_eProbe_temp_removed.snp'


if __name__ == '__main__':
    import argparse

    usage = '''python Snp_filter_bed.py -s SNPs_df.tsv -g chromosome_size.tsv -k keep_BED.tsv -d -50 ' \
            '-r remove_BED.tsv -m 50 -o output'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-s', '--snp_df', dest='snp_df', action='store', nargs='?',
                        help='Input SNPs data frame (tsv).',
                        metavar='FILE', required=True)
    parser.add_argument('-g', '--genome_size', dest='genome_size', action='store', nargs='?',
                        help='Genome size file (tsv).',
                        metavar='FILE', required=True)
    parser.add_argument('-k', '--k_bed', dest='k_bed', action='store', nargs='?',
                        help='Genome regions that need to retain their SNPs (bed format).',
                        metavar='FILE')
    parser.add_argument('-d', '--k_dist', dest='k_dist', action='store', nargs='?',
                        help='Distance around(+)/within(-) the regions needed to be kept (default: -50).',
                        default='-50', metavar='STRING')
    parser.add_argument('-r', '--r_bed', dest='r_bed', action='store', nargs='?',
                        help='Genome regions that need to retain their SNPs (bed format).',
                        metavar='FILE')
    parser.add_argument('-m', '--r_dist', dest='r_dist', action='store', nargs='?',
                        help='Distance around(+)/within(-) the regions needed to be removed (default: 50).',
                        default='50', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.',
                        metavar='STRING', required=True)
    args = parser.parse_args()

    temp_snp = SNP2BED(args.snp_df, args.output)

    if args.k_bed and check_bed(args.k_bed):
        print('Keeping SNPs in BED...', datetime.datetime.now())
        temp_bed = run_bedtools_slop(args.k_bed, args.genome_size, int(args.k_dist), args.output)
        kept_snp = run_bedtools_keep(temp_snp, temp_bed, args.output)
        print('Filtered out %s SNPs in this step.' % (diff_file_line_count(temp_snp, kept_snp)),
              datetime.datetime.now())
        mv_command = f'mv {kept_snp} {temp_snp}'
        subprocess.run(mv_command, shell=True, check=True)

    if args.r_bed and check_bed(args.r_bed):
        print('Removing SNPs in BED ...', datetime.datetime.now())
        temp_bed = run_bedtools_slop(args.r_bed, args.genome_size, int(args.r_dist), args.output)
        rm_snp = run_bedtools_remove(temp_snp, temp_bed, args.output)
        print('Filtered out %s SNPs in this step.' % (diff_file_line_count(temp_snp, rm_snp)),
              datetime.datetime.now())
        mv_command = f'mv {rm_snp} {temp_snp}'
        subprocess.run(mv_command, shell=True, check=True)

    awk_command = f'''awk 'BEGIN {{OFS="\\t"}} {{print $1, $2, $4, $5, $6}}' {temp_snp} > {args.output}.processing_SNPs.tsv'''
    subprocess.run(awk_command, shell=True, check=True)
    rm_command = f'rm {args.output}_eProbe_temp*'
    subprocess.run(rm_command, shell=True, check=True)
    before_filtering = str(int(subprocess.check_output(['wc', '-l', args.snp_df]).split()[0]))
    after_filtering = str(diff_file_line_count(args.snp_df, args.output + '.processing_SNPs.tsv'))
    print('Filtered out %s SNPs from %s SNPs in BED filtering.' % (after_filtering, before_filtering),
          datetime.datetime.now())
