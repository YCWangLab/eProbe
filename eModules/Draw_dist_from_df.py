import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def select_column_from_df(df, column_name):
    if column_name not in df.columns:
        raise ValueError(f'''Column '{column_name}' not found in the DataFrame.''')

    return df[column_name]


def plot_distribution(df_column, color, outpre, xlim, ylim, visual, style='darkgrid', stat='percent',
                      line_width=2, bins='auto',
                      label_fontsize=18, **kwargs):
    sns.set(style=style)

    plt.figure(figsize=(8, 6))

    line_kws = {'linewidth': line_width}

    line_kws.update(kwargs)

    if xlim != 'auto':
        bins = int((df_column.max() - df_column.min()) / float(bins))

    sns.histplot(df_column, alpha=.7, stat=stat, edgecolor=color, color=color, bins=bins, element=visual,
                 line_kws=line_kws)

    plt.title('Distribution Plot', fontsize=label_fontsize + 4, fontname='Arial', pad=20)
    plt.xlabel('Values', fontsize=label_fontsize + 2, fontname='Arial', labelpad=10)
    if stat == 'percent':
        plt.ylabel(f'{stat.capitalize()} (%)', fontsize=label_fontsize, fontname='Arial', labelpad=10)
    else:
        plt.ylabel(stat.capitalize(), fontsize=label_fontsize + 2, fontname='Arial', labelpad=10)

    if xlim != 'auto':
        plt.xlim([int(xlim.split(',')[0]), int(xlim.split(',')[1])])
    if ylim != 'auto':
        plt.ylim([float(ylim.split(',')[0]), float(ylim.split(',')[1])])

    plt.xticks(fontsize=label_fontsize, fontname='Arial')
    plt.yticks(fontsize=label_fontsize, fontname='Arial')

    plt.tight_layout()
    plt.savefig(f'{outpre}_{df_column.name}_dist.jpg', dpi=600)


if __name__ == '__main__':
    import argparse

    usage = '''python Draw_dist_from_df.py -d df.tsv -c col_name1 col_name2 -o prefix'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-d', '--df', dest='df', action='store', nargs='?',
                        help='Data frame.', metavar='FILE', required=True)
    parser.add_argument('-c', '--col_names', dest='col_name', action='store', nargs='+',
                        help='Column name(s) for plot distribution.', metavar='STRING', required=True)
    parser.add_argument('--sep', dest='separator', action='store', nargs='?', default='\t',
                        help='Separator of dataframe.', metavar='STRING')
    parser.add_argument('--colors', dest='colors', action='store', nargs='+',
                        help='Colors for plot (default: 5 default colors).',
                        default=None, metavar='STRING')
    parser.add_argument('--style', dest='style', choices=['darkgrid', 'whitegrid', 'dark', 'white', 'ticks'],
                        action='store', nargs='?', default='darkgrid',
                        help='Grid style [darkgrid, whitegrid, dark, white, ticks] for plot (default: darkgrid).',
                        metavar='STRING')
    parser.add_argument('--stat', dest='stat', action='store', nargs='?',
                        help='Stat style [percent, density, probability, count, frequency] for plot (default: percent).',
                        choices=['percent', 'density', 'probability', 'count', 'frequency'], default='percent',
                        metavar='STRING')
    parser.add_argument('--visual', dest='visual', action='store', nargs='?',
                        help='Visualization style (bar/step/poly, default: step).',
                        default='step', metavar='STRING')
    parser.add_argument('--bins', dest='bins', action='store', nargs='?',
                        help='Bin interval for histogram (default: auto).',
                        default='auto', metavar='STRING')
    parser.add_argument('--line_width', dest='line_width', action='store', nargs='?', default='2',
                        help='Line width for KDE plot (default: 2).', metavar='STRING')
    parser.add_argument('--xlim', dest='xlim', action='store', nargs='?',
                        help='Lower and upper limits of X axis, separated by comma (default: auto).',
                        default='auto', metavar='STRING')
    parser.add_argument('--ylim', dest='ylim', action='store', nargs='?',
                        help='Lower and upper limits of Y axis, separated by comma (default: auto).',
                        default='auto', metavar='STRING')
    parser.add_argument('--fontsize', dest='fontsize', action='store', nargs='?', default='18',
                        help='Fontsize (default: 18).', metavar='STRING')
    parser.add_argument('-o', '--output', dest='output', action='store', nargs='?',
                        help='Output prefix.', metavar='STRING', required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.df, sep=args.separator, header=0)
    if not args.colors:
        colors = ['#83639F', '#EA7827', '#C22f2F', '#449945', '#1F70A9']
    else:
        colors = args.colors

    for index, col in enumerate(args.col_name):
        df_col = select_column_from_df(df, col)

        if args.bins == 'auto':
            if col.lower() == 'gc':
                bins = 2
            elif col.lower() == 'tm':
                bins = 1
            elif col.lower() == 'complexity':
                bins = 0.1
            elif col.lower() == 'hairpin':
                bins = 2
            elif col.lower() == 'dimer':
                bins = 0.01
            else:
                bins = args.bins
        else:
            bins = args.bins

        if args.xlim == 'auto':
            if col.lower() in ['gc', 'tm']:
                xlim = '0,100'
            elif col.lower() == 'complexity':
                xlim = '0,5'
            elif col.lower() == 'hairpin':
                xlim = '0,100'
            elif col.lower() == 'dimer':
                xlim = '0,1'
            else:
                xlim = args.xlim
        else:
            xlim = args.xlim

        if args.ylim == 'auto':
            if col in ['gc', 'tm', 'complexity', 'hairpin']:
                if args.stat == 'percent':
                    ylim = '0,40'
                elif args.stat == 'density':
                    ylim = '0,0.4'
                else:
                    ylim = args.ylim
            elif col == 'dimer':
                if args.stat == 'percent':
                    ylim = '0,16'
                elif args.stat == 'density':
                    ylim = '0,0.16'
                else:
                    ylim = args.ylim
            else:
                ylim = args.ylim
        else:
            ylim = args.ylim

        color = colors[min(index, len(args.col_name))]
        plot_distribution(df_col, color, args.output, xlim, ylim, args.visual, style=args.style,
                          stat=args.stat,
                          line_width=args.line_width,
                          bins=bins, label_fontsize=int(args.fontsize))
