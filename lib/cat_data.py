# Script for untarring and concatenating the meta and sequence data
# downloaded from GISAID. 

import sys
import subprocess
import re
import os
import pandas as pd

in_dir = sys.argv[1]
metadata_list = []
sequence_list = []

def cat_raw(in_dir):
    metadata_list = []
    sequence_list = []
    for root,dirs,files in os.walk(in_dir):
        for name in files:
            if re.search(".tar", name):
                file = os.path.join(root, name)
                subprocess.run(["tar", "-xvf", file, "-C", root])

    for root,dirs,files in os.walk(in_dir):
        for name in files:
            if re.search(".tsv", name):
                file = os.path.join(root, name)
                metadata_list.append(file)
            if re.search(".fasta", name):
                file = os.path.join(root, name)
                with open(file, "r") as f:
                    sequence_list.append(f.readlines())

    all_meta = pd.read_csv(metadata_list[0], sep="\t")
    for meta_file in metadata_list[1:]:
        meta = pd.read_csv(meta_file, sep="\t")
        all_meta = pd.concat([all_meta, meta], ignore_index=True)

    os.makedirs("results/sequences")
    os.mkdir("results/metadata")

    with open("results/sequences/all_seq_concat.fasta", "w") as f:
        for seq in sequence_list:
            f.writelines(seq)

    all_meta.to_csv("results/metadata/all_meta_concat.csv")