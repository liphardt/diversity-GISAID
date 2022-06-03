source ~/anaconda3/etc/profile.d/conda.sh
conda activate nextstrain

mkdir aligned_fastas

for fasta in $(find week_subset -name "*fasta")
do
	name=$(echo $fasta | cut -d "/" -f 3 | cut -d "_" -f 1,2)
	echo $name
	nextalign --reference $1\
 	--sequences $fasta\
	--output-dir aligned_fastas\
	--output-basename ${name}
done
