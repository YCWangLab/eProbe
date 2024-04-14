import pandas as pd


def random_subsample_df(df, target_rows, seed):
    if not isinstance(target_rows, int) or target_rows <= 0:
        raise ValueError("Target rows for subsampling must be a positive integer")

    if not isinstance(seed, int):
        raise ValueError("Seed for sampling must be an integer")

    current_rows = df.shape[0]

    if current_rows <= target_rows:
        return df

    sample_ratio = target_rows / current_rows

    downsampled_df = df.sample(frac=sample_ratio, random_state=seed)

    return downsampled_df.sort_index()


def random_subsample_dfs(df1, df2, target_rows, seed):
    if not isinstance(target_rows, int) or target_rows < 0:
        raise ValueError("Target rows for subsampling must be a positive integer")
    if not isinstance(seed, int):
        raise ValueError("Seed for sampling must be an integer")

    inner_df = pd.merge(df1, df2, how='inner')

    df1_ = df1.loc[~df1.index.isin(inner_df.index)]
    df2_ = df2.loc[~df2.index.isin(inner_df.index)]

    total_rows = len(df1) + len(df2) - len(inner_df)  # Adjust total to avoid double-counting inner rows
    target_inner_rows = int(target_rows * (len(inner_df) / total_rows))
    target_df1_rows = int(target_rows * (len(df1_) / total_rows))
    target_df2_rows = target_rows - target_inner_rows - target_df1_rows  # Ensure total target_rows is met

    subsample_inner_df = inner_df.sample(n=target_inner_rows, random_state=seed)
    subsample_df1 = df1_.sample(n=target_df1_rows, random_state=seed)
    subsample_df2 = df2_.sample(n=target_df2_rows, random_state=seed)

    return pd.concat([subsample_df1, subsample_df2, subsample_inner_df], ignore_index=True).sort_index()


if __name__ == '__main__':
    import argparse

    usage = """python SNPs_subsampler.py -s 123 -r 20000 -o *.txt --df1 snp_df1 --df2 snp_df2 (optional)"""

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-s", "--seed", dest="seed", action="store", nargs='?', default="123",
                        help="Seed number for random subsampling.", metavar="STRING")
    parser.add_argument("-r", "--row", dest="rows", action="store", nargs='?', help="Target rows.",
                        metavar="STRING", required=True)
    parser.add_argument("--df1", dest="SNPs_df1", action="store",
                        help="SNPs dataframe1 for subsample.", metavar="FILE", required=True)
    parser.add_argument("--df2", dest="SNPs_df2", action="store", nargs='?',
                        help="SNPs dataframe2 for subsample (optional).", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)
    args = parser.parse_args()

    SNPs_df1 = pd.read_csv(args.SNPs_df1, sep='\t', header=0)
    if args.SNPs_df2:
        SNPs_df2 = pd.read_csv(args.SNPs_df2, sep='\t', header=0)
        out_df = random_subsample_dfs(SNPs_df1, SNPs_df2, int(args.rows), int(args.seed))
    else:
        out_df = random_subsample_df(SNPs_df1, int(args.rows), int(args.seed))

    out_df.to_csv(f"{args.output}.subsampled_{args.rows}.tsv", "\t", index=False)
