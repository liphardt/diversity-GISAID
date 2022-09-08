import os
import subprocess
import re
from Bio import AlignIO
from pathlib import Path


def align_seq(in_dir):
	os.makedirs("results/aligned_fastas")
	os.makedirs("results/other_next_out")
	out_fasta_dir = "results/aligned_fastas"
	out_other_dir = "results/other_next_out"
	root_path = Path(__file__).absolute()
	ref = root_path.parent / ".." / "reference/GCF_009858895.2_ASM985889v3_genomic.fna"
	for root,dirs,files in os.walk(in_dir):
		for names in files:
			if re.search(".fasta", names):
				full_path = os.path.join(root, names)
				name = '_'.join(names.split('_')[:2])
				subprocess.run(["nextalign", "--reference", ref, "--sequences",
				full_path, "--output-fasta", f"{out_fasta_dir}/{name}.fasta",
				"--output-insertions", f"{out_other_dir}/{name}.insertions.csv",
				"--output-errors", f"{out_other_dir}/{name}.errors.csv"])

def replace_codes(in_dir):
	iupac_ambig = ["R", "M", "W", "S", "Y", "K", "V", "H", "D", "B"]
	for root,dirs,files in os.walk(in_dir):
		for file in files:
			full_path = os.path.join(root,file)
			alignment = AlignIO.read(full_path, "fasta")
			for record in alignment:
				for letter in iupac_ambig:
					record.seq = record.seq.replace(letter, "N")
			AlignIO.write(alignment, full_path, "fasta")
