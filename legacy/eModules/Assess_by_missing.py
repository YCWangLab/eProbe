import random
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysam
import itertools
from pipe_controler import delete_files


def simulate_random_missing(vcf_file, missing_rates, missing_genotype, output_file):
    with open(vcf_file, 'r') as file:
        with open(output_file, 'w') as outfile:
            for line in file:
                if line.startswith('#'):
                    outfile.write(line)
                else:
                    variant_data = line.strip().split('\t')
                    info_columns = variant_data[:9]
                    genotype_columns = variant_data[9:]

                    simulated_genotypes = []
                    for genotype in genotype_columns:
                        if random.random() < random.choice(missing_rates):
                            simulated_genotypes.append(missing_genotype)
                        else:
                            simulated_genotypes.append(genotype)

                    simulated_line = '\t'.join(info_columns + simulated_genotypes) + '\n'
                    outfile.write(simulated_line)


def get_pca(vcf, output):
    plink_command = f'plink --allow-extra-chr --vcf {vcf} --make-bed --out {output} && ' \
                    f'plink --allow-extra-chr --bfile {output} --pca --out {output}_pca'
    subprocess.run(plink_command, shell=True, check=True)
    delete_files(
        [f'{output}.bed', f'{output}.bim', f'{output}.fam', f'{output}.log', f'{output}.nosex',
         f'{output}_pca.nosex', f'{output}_pca.log'])
    return f'{output}_pca.eigenvec', f'{output}_pca.eigenval'


def read_sample_ids(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    sample_ids = list(vcf.header.samples)
    vcf.close()
    return sample_ids


def read_pca(eigenvec, eigenval, new_ids):
    eigenvectors = pd.read_csv(eigenvec, delim_whitespace=True, header=None)
    eigenvectors.drop(columns=[1], inplace=True)
    num_columns = eigenvectors.shape[1]
    column_names = ['ID'] + [f'PC{i}' for i in range(1, num_columns)]
    eigenvectors.columns = column_names
    eigenvectors['ID'] = new_ids

    eigenvalues = pd.read_csv(eigenval, header=None)
    total_variance = eigenvalues[0].sum()
    contributions = (eigenvalues[0] / total_variance) * 100
    contributions_dict = {f'PC{i + 1}': contribution for i, contribution in enumerate(contributions)}

    return eigenvectors, contributions_dict


def read_pop(pop_file):
    pop_data = pd.read_csv(pop_file, delimiter='\t', header=None, names=['ID', 'Population'])
    return pop_data


def pca_draw2(eigenvectors, contributions_dict, x, y, output, pop_file=None):
    fig, ax = plt.subplots(figsize=(8, 6))
    if pop_file is not None:
        pop_data = read_pop(pop_file)
        if len(pop_data) != len(eigenvectors):
            print('The number of samples in VCF does not match the number in the population file. '
                  'Only the intersect of these two will be used!')
        merged_df = pd.merge(eigenvectors, pop_data, on='ID', how='inner')
        populations = merged_df['Population'].unique()
        color_map = plt.get_cmap('Dark2')
        colors = color_map(np.linspace(0, 1, len(populations)))
        population_colors = {pop: colors[i] for i, pop in enumerate(populations)}

        for pop in populations:
            subset = merged_df[merged_df['Population'] == pop]
            ax.scatter(subset[f'PC{int(x)}'], subset[f'PC{int(y)}'], color=population_colors[pop], label=pop,
                       alpha=0.8, s=30, edgecolors='black', linewidths=0.7)

        ax.legend(title='Population', loc='best', fontsize=10, frameon=False)

    else:
        num_points = eigenvectors.shape[0]
        np.random.seed(1)
        colors = np.random.rand(num_points)
        ax.scatter(eigenvectors[f'PC{int(x)}'], eigenvectors[f'PC{int(y)}'], c=colors, cmap='Dark2', alpha=0.8, s=30,
                   edgecolors='black', linewidths=0.7)

    ax.set_xlabel(f'''PC{int(x)} {round(contributions_dict[f'PC{int(x)}'], 2)}%''', fontsize=18, fontname='Arial')
    ax.set_ylabel(f'''PC{int(y)} {round(contributions_dict[f'PC{int(y)}'], 2)}%''', fontsize=18, fontname='Arial')
    ax.tick_params(axis='both', which='major', labelsize=15)
    plt.tight_layout()
    plt.savefig(f'''{output}_{f'PC{int(x)}'}_{f'PC{int(y)}'}.jpg''', dpi=800)


if __name__ == '__main__':
    import argparse

    usage = '''python Assess_by_missing.py -v *.vcf -r 0.2,0.4,0.6 -g ./.:0,0,0 -o simulated'''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-v', '--vcf', dest='vcf', action='store', nargs='?',
                        help='VCF file (uncompressed)',
                        metavar='FILE', required=True)
    parser.add_argument('-r', '--missing_rate', dest='missing_rate', action='store', nargs='?',
                        help='Missing rate. Can be fixed for all individuals (e.g., 0.2), '
                             'or a list for random selection (default: 0.2,0.4,0.6)',
                        default='0.2,0.4,0.6', metavar='STRING')
    parser.add_argument('-g', '--missing_genotype', dest='missing_genotype', action='store', nargs='?',
                        help='''Genotype for replacement, normally './.:0,0,0' for VCF called from BCFtools'
                             'and '.:0,0:0:.:0,0' for GATK)''',
                        default='./.:0,0,0', metavar='STRING')
    parser.add_argument('--pc', dest='pc', action='store', nargs='?',
                        help='How many PCs to draw (default: 1,2).',
                        default='1,2', metavar='FILE')
    parser.add_argument('--pop_data', dest='pop_data', action='store', nargs='?',
                        help='A population information to color the map.'
                             'Should be a tab-separated file and each line consist of '
                             'Sample_ID and Population_ID, without header.',
                        default=None, metavar='FILE')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.', metavar='STRING', required=True)
    args = parser.parse_args()

    missing_rate = list(map(float, args.missing_rate.split(',')))
    simulated_vcf = f'{args.output}.simulated.vcf'
    simulate_random_missing(args.vcf, missing_rate, args.missing_genotype, simulated_vcf)

    origin_vec, origin_val = get_pca(args.vcf, args.output + '_original')
    simulated_vec, simulated_val = get_pca(simulated_vcf, args.output + '_simulated')

    origin_vec_df, origin_val_dict = read_pca(origin_vec, origin_val, read_sample_ids(args.vcf))
    simulated_vec_df, simulated_val_dict = read_pca(simulated_vec, simulated_val, read_sample_ids(simulated_vcf))

    PCs = list(itertools.combinations(args.pc.split(','), 2))

    for PC in PCs:
        pca_draw2(origin_vec_df, origin_val_dict, PC[0], PC[1], args.output + '_original', pop_file=args.pop_data)
        pca_draw2(simulated_vec_df, simulated_val_dict, PC[0], PC[1], args.output + '_simulated',
                  pop_file=args.pop_data)
