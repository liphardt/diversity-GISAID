from typing import SupportsRound
import lib.cat_data as cd
import lib.preprocess as prepro
import lib.align as align
import lib.group_cases as gc
import argparse
import pandas as pd
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--inputdir', help = 'Raw input directory from GISAID')
parser.add_argument('-S', '--state', help = 'Which state you are processing. For subsetting case data.')
args = parser.parse_args()

root_dir = Path(__file__).absolute()
r_script = f"{root_dir.parent}/lib/plot_covid.R"

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
    align.align_seq("results/week_subset")
    print("Replacing IUPAC ambiguity codes")
    align.replace_codes("results/aligned_fastas")
    print(f"Collecting case data for {args.state}")
    print("To keep case data updated go to")
    print("https://data.cdc.gov/Case-Surveillance/COVID-19-Case-Surveillance-Public-Use-Data/vbim-akqf/data")
    print("to download latest case data.")
    gc.group_cases_by_state(args.state)
    print("Finished grouping state data.")
    subprocess.run("sed -i '/^[^>]/ s/-/N/g' results/aligned_fastas/*",shell=True)
    print("Calculating variant replacement and nucleotide diversity through time.")
    r_args = ["Rscript", r_script, f"results/{args.state}_processed_metadata.csv", f"results/{args.state}_all_case_data.csv", f"results/{args.state}_case_period_counts.csv", args.state]
    print(r_args)
    subprocess.run(r_args)
