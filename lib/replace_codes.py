import sys
from Bio import AlignIO
import re

alignment = AlignIO.read(sys.argv[1], "fasta")

iupac_ambig = ["R", "M", "W", "S", "Y", "K", "V", "H", "D", "B"]


print("Such a pain in the ass...")
for record in alignment:
	for let in iupac_ambig:
		if let in record.seq:
			record.seq = record.seq.replace(let, "N")
print("Okay, writing the new one...")
AlignIO.write(alignment, sys.argv[1], "fasta")
