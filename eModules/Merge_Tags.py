import pandas as pd


def merge_2_df(df1, df2, how):
    merged_df = pd.merge(df1, df2, on=list(set(df1.columns) & set(df2.columns)), how=how)
    return merged_df


def merge_dfs(df_list, how):
    merged_df = merge_2_df(df_list[0], df_list[1], how)
    for df in df_list[2:]:
        merged_df = merge_2_df(merged_df, df, how)

    return merged_df


def read_txt_files_to_df(file_paths):
    df_list = []
    for file_path in file_paths:
        try:
            df = pd.read_csv(file_path, sep='\t', header=0)
            df_list.append(df)
        except Exception as e:
            print(f"Reading {file_path} failed: {e}")

    return df_list


if __name__ == '__main__':
    import argparse

    usage = """ python Merge_Tags.py --how outer -o output -f df1 df2 ..."""

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("--how", dest="how", action="store", nargs='?', default="outer",
                        help="The method [inner, outer, left, right] to merge two data frames (default: outer).",
                        metavar="STRING")
    parser.add_argument("-f", "--files", dest="files", action="store", nargs='+',
                        help="Tsv files to merge.", metavar="FILE")
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)
    args = parser.parse_args()

    df_list = read_txt_files_to_df(args.files)

    merged_df = merge_dfs(df_list, args.how)

    merged_df.to_csv(args.output + "_merged_tags.tsv", sep='\t', index=False)
