#!/bin/bash --login

######################
# Author: Benjamin Anderson
# Date: June 2024
# Description: translate multifastas that are in-frame, then align with MAFFT and convert back to nucleotide  
# Note: requires a Singularity container (phylo.sif) with the dependencies (MAFFT and pal2nal)
#	It requires two Python scripts for removing entries (optional) and for translating
#	It assumes the input multifastas are named *_name.fasta (so, will keep the last field separated by "_"
#	or everything if no underscore)
#	Arguments are -f for the input files directory (will run on all fasta files there)
#	and -d for the dropfile of samples to exclude (one per line; optional)
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=transalign
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=transalign_%j.log


# set the container and scripts
container="/software/projects/pawsey0220/anderson/singularity-containers/phylo.sif"
remove_script="/software/projects/pawsey0220/anderson/scripts/remove_fastas.py"
translate_script="/software/projects/pawsey0220/anderson/scripts/translate.py"


# load modules
module load python/3.11.6 singularity/4.1.0-slurm


# parse the command line
while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-d)
	dropfile="$(readlink -f $2)"
	shift
	shift
	;;
	-f)
	fasta_dir="$(readlink -f $2)"
	shift
	shift
	;;
esac
done

# check args
if [ -z "$dropfile" ]; then
	drop="no"
else
	drop="yes"
fi

if [ -z "$fasta_dir" ]; then
	echo -e "\nPlease specify where the fasta files are\n"
	exit 1
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting translation and alignment at $(date)\n**********"


# if samples need to be dropped, do so, otherwise copy files to current directory
if [ "$drop" = "yes" ]; then
	for fasta_file in "$fasta_dir"/*.fasta; do
		python3 "$remove_script" -f "$dropfile" "$fasta_file"
		for file in mod*.fasta; do 
			newname=$(echo $file | awk -F_ '{print $NF}'); mv $file $newname
		done
	done
else
	for fasta_file in "$fasta_dir"/*.fasta; do
		newname=$(basename $fasta_file | awk -F_ '{print $NF}')
		cp "$fasta_file" ./"$newname"
	done
fi


# translate the files and align them
python3 "$translate_script" *.fasta
sed -i 's/\*/X/g' *_prot.fasta		# convert stop codons to "X"
for translation in *_prot.fasta; do
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
		mafft --auto --thread "$SLURM_CPUS_PER_TASK" "$translation" > "${translation/\.fasta/\.aln}"
done


# convert back to nucleotide
for alignment in *_prot.aln; do
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec "$container" \
		pal2nal "$alignment" "${alignment/_prot\.aln/\.fasta}" \
		-output fasta > "${alignment/prot\.aln/exon_aligned\.fasta}"
done


# rename final files to locus names and remove intermediate files
for file in *exon_aligned.fasta; do
	mv "$file" "${file/_exon_aligned/}"
	rm "${file/exon_aligned\.fasta/prot}"*
done


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished translation and alignment run at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours"
