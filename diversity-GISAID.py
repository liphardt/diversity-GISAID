import lib.cat_data as cd
import lib.preprocess as prepro
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--inputdir', help = 'Raw input directory from GISAID')
parser.add_argument('-S', '--state', help = 'Which state you are processing. For subsetting case data.')
args = parser.parse_args()


if __name__ == "__main__":
    print(f"Untarring and concatenating GISAID data in {args.inputdir}")
    cd.cat_raw(args.inputdir)
    print("Finished concatenating raw data")
    print("Processing GISAID data...")
    gis_meta = pd.read_csv("results/metadata/all_meta_concat.csv")
    gis_sequences = "results/sequences/all_seq_concat.fasta"
    fixed_meta = prepro.fix_time(gis_meta)
    prepro.add_data_columns(fixed_meta)
    prepro.count_total(fixed_meta)
    prepro.check_loc(fixed_meta)
    prepro.add_who(fixed_meta)
    outmeta = f"results/{args.state}_processed_metadata.csv"
    print(f"Writing modified metadata to file {outmeta}")
    fixed_meta.to_csv(outmeta, index = False)
    print("Pulling metadata and sequences for week intervals starting at 1 January, 2020")
    prepro.pull_meta_seqs(fixed_meta, gis_sequences)
    print("Finished data preprocessing")
    print("Aligning sequences to SARS reference")
    