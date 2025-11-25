#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: June 2024
# Modified: May 2025 (upgrade HybPiper v. 2.3.2, with changes to output)
# Description: combine results from multiple HybPiper runs  
# Note: submit this from the directory where the folder "results" is
#	It requires two Python scripts for combining files
#	If you want to delete the individual files, set argument `-d "yes"`
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=combHP
#SBATCH --partition=copy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=combHP_%j.log


# load the module
module load python/3.11.6


# set the scripts
len_script="/software/projects/pawsey0220/anderson/scripts/hybpiper_loci_lengths.py"
comb_script="/software/projects/pawsey0220/anderson/scripts/fasta_combining.py"


# parse the command line
while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-d)
	delete="$2"
	shift
	shift
	;;
esac
done

# check arg
if [ -z "$delete" ]; then
	delete="no"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nCombining results at $(date)\n**********"


# change into the results directory
cd results


# combine stats files into one
ls -d */ | sed 's/\///' > temp_samples.txt
sample=$(head -n 1 temp_samples.txt)
head -n 1 "$sample"/*stats.tsv > combined_stats.tsv
for samp in $(cat temp_samples.txt); do 
	tail -1 "$samp"/*stats.tsv >> combined_stats.tsv
done


# combine lengths files into one
# NOTE: the gene lengths are in a different order for each sample
# NOTE: the mean lengths are incorrect (different for each)
python3 "$len_script" */*lengths.tsv


# combine the fasta files into one per locus
if [ ! -d dna_seqs ]; then
	mkdir dna_seqs
fi
if [ ! -d paralog_seqs ]; then
	mkdir paralog_seqs
fi
if [ ! -d supercontig_seqs ]; then
	mkdir supercontig_seqs
fi


# create a function for combining files
combine_fastas () {
	loc="$1"
	python3 "$comb_script" -o "$loc"_temp.fasta */dna_seqs/"$loc".FNA
	mv "$loc"_temp.fasta dna_seqs/"$loc".fasta

	python3 "$comb_script" -o "$loc"_temp2.fasta */paralogs_all/"$loc"_paralogs_all.fasta
	mv "$loc"_temp2.fasta paralog_seqs/"$loc".fasta

	python3 "$comb_script" -o "$loc"_temp3.fasta */supercontig_seqs/"$loc".fasta
	mv "$loc"_temp3.fasta supercontig_seqs/"$loc".fasta
}


# to parallelise, split the loci into groups
ls */dna_seqs/*.FNA | cut -f 3 -d "/" | cut -f 1 -d "." | sort | uniq > temp_loci.txt
split -l "$SLURM_CPUS_PER_TASK" temp_loci.txt temp_subset_
for subfile in temp_subset_*
do
	for loc in $(cat $subfile)
	do
		combine_fastas "$loc" &
	done
	wait
done


# remove the individual sequences
if [ $delete == "yes" ]; then
	find -P */dna_seqs/ -type f -print0 -o -type l -print0 | xargs -0 munlink
	find -P */paralogs_all/ -type f -print0 -o -type l -print0 | xargs -0 munlink
	find -P */paralogs_no_chimeras/ -type f -print0 -o -type l -print0 | xargs -0 munlink
	find -P */supercontig_seqs/ -type f -print0 -o -type l -print0 | xargs -0 munlink

	find -P */dna_seqs/ -type d -empty -delete
	find -P */paralogs_all/ -type d -empty -delete
	find -P */paralogs_no_chimeras/ -type d -empty -delete
	find -P */supercontig_seqs/ -type d -empty -delete
fi
rm temp*


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished combining results after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours"
