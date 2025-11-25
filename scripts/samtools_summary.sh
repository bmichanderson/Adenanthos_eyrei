#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Modified: Mar 2025 (removed module load; adjusted -e to allow specifying number of bp extended; added output header)
# Description: run a Samtools summary job for mapped reads on Setonix
# Note: this will report read depth as well as output a consensus including ambiguities
#	It should be called by the array script `samtools_array.sbatch` for a group of samples
#
#	For the purposes of stats (length, mean depth, etc.), het counts and the consensus sequence,
#	this script can exclude the first and last specified bp of sequence (turn on by providing -e {length} [default 0])
#	The only other argument is for the name of the sample (-s), which should match a folder in the working directory
#
#	If interested, Samtools arguments for consensus generation can be adjusted manually:
#	minimum read depth = how many of the top base + second most common base must be present to make a call [default 15],
#	het fraction = proportion of top base needed to have a second base for reporting ambiguities [default 0.4], and
#	call fraction = minimum proportion of bases used in the call to make a call [default 0.6]
#
#	As I understand it, given a position has 15 read depth (9 T, 4 A, 2 C), the call would be "N"
#	The second highest base is included for ambiguities because 4/9 > 0.4 (ambiguity would be "W")
#	The call is "N" because top + second base (9 + 4) is < 15; if min depth is lowered to 13, then the call is "W"
#	Call fraction is satisfied as well (9 + 4) / 15 > 0.6
#	The weird case would be (I think) if the bases were 9 T, 3 A, 3 C
#	het fraction 3/9 = 0.33 < 0.4, so no second base; as a result, depth is 9, which is below min depth and call = "N", not "T"
#
#	This script will create a couple of files, one of which is called `<sample>_samtools_summary.txt`
#	and contains columns:
#		sample	length	min_depth	max_depth	mean_depth stdev_depth	count_hets	pos_hets	count_Ns	pos_Ns
######################


# define container
container="/software/projects/pawsey0220/anderson/singularity-containers/samtools.sif"


# define Samtools parameters
min_depth="15"
het_fraction="0.4"
call_fraction="0.6"


# Load the Singularity module (may need to update depending on Setonix updates)
module load singularity/4.1.0-slurm


# parse the command line
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
	-e)
	extend="$2"
	shift
	shift
	;;
	-s)
	sample="$2"
	shift
	shift
	;;
esac
done


# check args
if [ -z "$sample" ]; then
	echo -e "\nPlease specify a sample\n"
	exit 1
fi
if [ -z "$extend" ]; then
	extend=0
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting Samtools summary at $(date)\n**********"


# Run the depth summary and consensus generation for each sample
if [ -d "$sample" ]; then
	cd "$sample"
	echo -e "\nRunning depth assessment for sample $sample\n"
	singularity exec -H "$(pwd)" "$container" samtools depth -a *.bam > "$sample"_depth.txt
	if [ "$extend" -gt 0 ]; then
		# to calculate depth stats, exclude specified bp from each end of the mapping to account for drop off
		cut -f 3 "$sample"_depth.txt | tail -n +$((extend + 1)) | head -n -"$extend" > temp.txt
		length=$(cat temp.txt | wc -l)
		minimum=$(sort -n temp.txt | head -n 1)
		maximum=$(sort -n -r temp.txt | head -n 1)
		mean=$(awk '{ sum += $1 } END { print sum / NR }' temp.txt)
		stdev=$(awk -v m="$mean" '{ calc = ($1 - m)^2; sum += calc } END { print sqrt( sum / NR ) }' temp.txt)

		echo -e "\nRunning consensus generation for sample $sample\n"
		singularity exec -H "$(pwd)" "$container" samtools consensus --mode simple -a \
			--ambig --min-depth "$min_depth" --het-fract "$het_fraction" --call-fract "$call_fraction" \
			-o "$sample"_consensus.fasta *.bam
		sequence=$(grep -v ">" "$sample"_consensus.fasta | tr -d '\n')

		count_hets=$(echo "${sequence:extend:-extend}" | awk '{ print gsub(/[KMRSWYVHDB]/, "") }')
		echo "${sequence:extend:-extend}" | awk -v extend="$extend" '{ while (match($0, /[KMRSWYVHDB]/)) {
				var1 = RSTART + extend; var2 = substr($0, RSTART, 1); i = var2":"var1; arr[i]++;
				$0 = substr($0, 1, RSTART-1) " " substr($0, RSTART+1)
				}
			}
			{ for (i in arr) {printf "%s ", i}}' > temp.txt
		pos_hets=$(cat temp.txt)

		count_ns=$(echo "${sequence:extend:-extend}" | awk '{ print gsub(/N/, "") }')
		echo "${sequence:extend:-extend}" | awk -v extend="$extend" '{ while (match($0, /N/)) {
				var1 = RSTART + extend; var2 = substr($0, RSTART, 1); i = var2":"var1; arr[i]++;
				$0 = substr($0, 1, RSTART-1) " " substr($0, RSTART+1)
				}
			}
			{ for (i in arr) {printf "%s ", i}}' > temp.txt
		pos_ns=$(cat temp.txt)

		echo -e "sample\tlength\tmin_depth\tmax_depth\tmean_depth\tstdev_depth\tcount_hets\tpos_hets\tcount_Ns\tpos_Ns" \
			> "$sample"_samtools_summary.txt
		paste <(echo $sample) <(echo $length) <(echo $minimum) <(echo $maximum) <(echo $mean) \
			<(echo $stdev) <(echo $count_hets) <(echo $pos_hets) <(echo $count_ns) <(echo $pos_ns) \
			>> "$sample"_samtools_summary.txt
		rm temp.txt
		cd ..
	else
		cut -f 3 "$sample"_depth.txt > temp.txt
		length=$(cat temp.txt | wc -l)
		minimum=$(sort -n temp.txt | head -n 1)
		maximum=$(sort -n -r temp.txt | head -n 1)
		mean=$(awk '{ sum += $1 } END { print sum / NR }' temp.txt)
		stdev=$(awk -v m="$mean" '{ calc = ($1 - m)^2; sum += calc } END { print sqrt( sum / NR ) }' temp.txt)

		echo -e "\nRunning consensus generation for sample $sample\n"
		singularity exec -H "$(pwd)" "$container" samtools consensus --mode simple -a \
			--ambig --min-depth "$min_depth" --het-fract "$het_fraction" --call-fract "$call_fraction" \
			-o "$sample"_consensus.fasta *.bam
		sequence=$(grep -v ">" "$sample"_consensus.fasta | tr -d '\n')

		count_hets=$(echo "${sequence}" | awk '{ print gsub(/[KMRSWYVHDB]/, "") }')
		echo "${sequence}" | awk '{ while (match($0, /[KMRSWYVHDB]/)) {
				var1 = RSTART; var2 = substr($0, RSTART, 1); i = var2":"var1; arr[i]++;
				$0 = substr($0, 1, RSTART-1) " " substr($0, RSTART+1)
				}
			}
			{ for (i in arr) {printf "%s ", i}}' > temp.txt
		pos_hets=$(cat temp.txt)

		count_ns=$(echo "${sequence}" | awk '{ print gsub(/N/, "") }')
		echo "${sequence}" | awk '{ while (match($0, /N/)) {
				var1 = RSTART; var2 = substr($0, RSTART, 1); i = var2":"var1; arr[i]++;
				$0 = substr($0, 1, RSTART-1) " " substr($0, RSTART+1)
				}
			}
			{ for (i in arr) {printf "%s ", i}}' > temp.txt
		pos_ns=$(cat temp.txt)

		echo -e "sample\tlength\tmin_depth\tmax_depth\tmean_depth\tstdev_depth\tcount_hets\tpos_hets\tcount_Ns\tpos_Ns" \
			> "$sample"_samtools_summary.txt
		paste <(echo $sample) <(echo $length) <(echo $minimum) <(echo $maximum) <(echo $mean) \
			<(echo $stdev) <(echo $count_hets) <(echo $pos_hets) <(echo $count_ns) <(echo $pos_ns) \
			>> "$sample"_samtools_summary.txt
		rm temp.txt
		cd ..
	fi
else
	echo -e "\nSample $sample doesn't have a directory in the launch directory\n"
	exit 1
fi


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
echo -e "\nFinished Samtools summary at $(date) after running for" \
	"$duration seconds or $duration_mins minutes\n"
echo "Job ${SLURM_ARRAY_TASK_ID} for sample "$sample" finished in $duration seconds"
