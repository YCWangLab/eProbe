import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
import fastcluster
import subprocess
from pipe_controler import delete_files


def get_ibs(vcf, output):
    plink_command = f'plink --vcf {vcf} --make-bed --out {output} --allow-extra-chr --const-fid && ' \
                    f'plink --bfile {output} --allow-extra-chr --distance square 1-ibs --out {output}'
    subprocess.run(plink_command, shell=True, check=True)
    paste_command = f'paste {output}.mdist.id {output}.mdist | cut -f2- > {output}.mat'
    subprocess.run(paste_command, shell=True, check=True)
    delete_files(
        [f'{output}.bed', f'{output}.bim', f'{output}.fam', f'{output}.log', f'{output}.nosex',
         f'{output}.mdist', f'{output}.mdist.id'])
    return f'{output}.mat'


def plot_heatmap_with_clustering(matrix, filename, linewidths, cluster, cmap, vmax):
    if cluster == 'on':
        linkage = fastcluster.linkage_vector(matrix, method='ward')
        g = sns.clustermap(matrix, cmap=cmap, vmin=0, vmax=vmax, linewidths=linewidths,
                           figsize=(10, 10), cbar_kws={'shrink': .5}, annot=False,
                           xticklabels=False, yticklabels=False, row_linkage=linkage, col_linkage=linkage)
        g.ax_heatmap.set_ylabel('')
        g.savefig(filename, dpi=800)
        plt.close(g.fig)
    else:
        g = sns.heatmap(matrix, cmap=cmap, vmin=0, vmax=vmax, linewidths=linewidths,
                    annot=False, xticklabels=False, yticklabels=False, cbar=True)
        g.ax_heatmap.set_ylabel('')
        g.savefig(filename, dpi=800)
        plt.close()


def calculate_similarity_and_plot_heatmap(matrix1, matrix2, prefix, cmap, cluster=False, linewidths=0):
    max_value = max(matrix1.max().max(), matrix2.max().max())
    plot_heatmap_with_clustering(matrix1, f'{prefix}_All_distance.jpg', linewidths, cluster, cmap, vmax=max_value)
    plot_heatmap_with_clustering(matrix2, f'{prefix}_Part_distance.jpg', linewidths, cluster, cmap, vmax=max_value)


def cal_correlation(matrix1, matrix2, method):
    if method == 'Pearson':
        correlation = np.corrcoef(matrix1.values.flatten(), matrix2.values.flatten())[0, 1]
    elif method == 'Spearman':
        correlation = stats.spearmanr(matrix1.values.flatten(), matrix2.values.flatten())[0]
    else:
        raise ValueError(f'Unsupported correlation method "{method}". Supported methods are "Pearson" and "Spearman".')

    return correlation


def manhattan_distance(matrix1, matrix2):
    distance = np.sum(np.abs(matrix1.values.flatten() - matrix2.values.flatten()))
    return distance


if __name__ == '__main__':
    import argparse

    usage = '''python Assess_by_distance.py --vcf1 vcf1 --vcf2 vcf2 -o out_prefix'''
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('--vcf1', dest='vcf1', action='store', nargs='?',
                        help='Complete vcf file.', metavar='FILE')
    parser.add_argument('--vcf2', dest='vcf2', action='store', nargs='?',
                        help='Subset vcf file.', metavar='FILE')
    parser.add_argument('--cluster', dest='cluster', action='store', nargs='?',
                        help='Conduct hierarchical clustering (on/off, default: on).',
                        default='on', metavar='STRING')
    parser.add_argument('--cor', dest='cor', action='store', nargs='?',
                        help='Method for calculating correlation coefficient (default: Pearson).',
                        choices=['Pearson', 'Spearman'], default='Pearson', metavar='STRING')
    parser.add_argument('--cmap', dest='cmap', action='store', nargs='?',
                        help='Cmap [try flare, crest, viridis_r] for plotting heatmap (default: YlOrBr). ',
                        default='YlOrBr', metavar='STRING')
    parser.add_argument('--index_row', dest='index_row', action='store', nargs='?',
                        help='Whether contain header row (True/False, default: False).',
                        choices=[True, False], default=False, metavar='STRING')
    parser.add_argument('--index_col', dest='index_col', action='store', nargs='?',
                        help='Whether contain index column (True/False, default: True).',
                        default=True, choices=[True, False], metavar='STRING')
    parser.add_argument('-o', '--out', dest='output', action='store', nargs='?',
                        help='Prefix for output.', metavar='STRING')
    args = parser.parse_args()

    print('Calculating IBS distance ...')
    matrix1 = get_ibs(args.vcf1, f'{args.output}_All')
    matrix2 = get_ibs(args.vcf2, f'{args.output}_Part')

    if args.index_row:
        header = 0
    else:
        header = None

    if args.index_col:
        index_col = 0
    else:
        index_col = None

    matrix1 = pd.read_csv(matrix1, delimiter='\t', index_col=index_col, header=header)
    matrix2 = pd.read_csv(matrix2, delimiter='\t', index_col=index_col, header=header)

    correlation = cal_correlation(matrix1, matrix2, args.cor)
    print(f'{args.cor} correlation: {str(correlation)}')

    manhattan_dis = manhattan_distance(matrix1, matrix2)
    print(f'Manhattan distance: {str(manhattan_dis)}')

    print('Plotting heatmap ...')
    calculate_similarity_and_plot_heatmap(matrix1, matrix2, args.output, args.cmap, cluster=args.cluster)
