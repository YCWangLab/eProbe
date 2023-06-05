import random
import argparse
usage = """python Assess_By_Missing.py -v *.vcf -r 0.2,0.4,0.6 -g ./.:0,0,0 -o *.simulate.vcf"""
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-v", "--vcf", dest="vcf", action="store", nargs='?',
                    help="VCF file", metavar="FILE")
parser.add_argument("-r", "--missing_rate", dest="missing_rate", action="store", nargs='?',
                    help="Missing rate. Can be fixed for all individuals (e.g., 0.2), or a list for random selection (e.g., 0.2,0.4,0.6)",
                    metavar="FILE")
parser.add_argument("-g", "--missing_genotype", dest="missing_genotype", action="store", nargs='?',
                    help="Genotype for replacement, normally './.:0,0,0' for VCF called from BCFtools and '.:0,0:0:.:0,0' for GATK)", metavar="FILE")
parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                    help="Output file name", metavar="FILE")
args = parser.parse_args()
vcf_file = args.vcf
miss_r = list(map(float, args.missing_rate.split(",")))
miss_g = args.missing_genotype
output = args.output

def simulate_random_missing(vcf_file, missing_rate, missing_genotype):
    # reading vcf
    with open(vcf_file, 'r') as file:
        lines = file.readlines()

    # parsing header
    header_lines = []
    variant_lines = []
    for line in lines:
        if line.startswith('#'):
            header_lines.append(line)
        else:
            variant_lines.append(line)

    # parsing sample names
    sample_names = []
    for line in header_lines:
        if line.startswith('#CHROM'):
            sample_names = line.strip().split('\t')[9:]
            break

    # missing rate for each sample
    missing_counts = []
    missing_rates = []
    if len(missing_rate) == 1:
        for sample in sample_names:
            missing_count = int(len(variant_lines) * missing_rate[0])
            missing_counts.append(missing_count)
            missing_rates.append(missing_rate[0])
    elif len(missing_rate) >= 1:
        for sample in sample_names:
            indv_missing_rate = random.choice(missing_rate)
            missing_rates.append(indv_missing_rate)
            missing_count = int(len(variant_lines) * indv_missing_rate)
            missing_counts.append(missing_count)

    # missing simulation for each SNPs
    simulated_variants = []
    for line in variant_lines:
        variant_data = line.strip().split('\t')
        simulated_data = variant_data[:9]  # the first 9 columns of information
        # missing simulation for each sample in a SNPs
        for i, sample in enumerate(sample_names):
            missing_count = missing_counts[i]
            genotype = variant_data[i + 9]
            # random missing based on missing rate
            if missing_count > 0 and random.random() < missing_rates[i]:
                genotype = missing_genotype
                missing_count -= 1

            simulated_data.append(genotype)
            missing_counts[i] = missing_count

        simulated_variants.append('\t'.join(simulated_data))

    # construct simulated vcf content
    simulated_vcf_content = ''.join(header_lines) + '\n' + '\n'.join(simulated_variants)

    # writing
    with open(output, 'w') as file:
        file.write(simulated_vcf_content)

