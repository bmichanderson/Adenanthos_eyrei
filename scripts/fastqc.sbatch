#!/bin/bash -l

######################
# Author: B. Anderson
# Date: Mar 2022
# Modified: June 2022, Oct 2023, Aug 2024
# Description: run FastQC analyses on files
# Note: this will create a "fastqc" folder in the current directory if one does not exist
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=fastqc
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=10
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --output=fastqc_%j.log


# a function for when the script is called incorrectly or without arguments
usage()
{
	echo -e "This script is for running FastQC analyses on Setonix\n"
	echo -e "Usage:\tsbatch fastqc.sbatch file1 file2 ...\n"
	exit 1
}


# load the module
module load fastqc/0.11.9--hdfd78af_1


# determine how many files to process at the same time
num_per_batch=$((SLURM_CPUS_PER_TASK / 2))


if [ $# -eq 0 ]; then		# if there are no command arguments
	usage
fi


# Record and echo start time
start="$(date +%s)"
echo -e "\nStarting FastQC analysis at $(date)\n"


# check if the output directory exists and if not, create it
if [ ! -d fastqc ]; then
	mkdir fastqc
fi


# run FastQC on your input files
echo $* | xargs -n "$num_per_batch" > temp_job_lines.txt

num_jobs=$(cat temp_job_lines.txt | wc -l)
echo -e "Will run $num_jobs FastQC analyses\n"
max_jobs="$SLURM_NTASKS"
curr_jobs=0
for ((job=1; job<=(($num_jobs)); job++)); do
	((curr_jobs >= max_jobs)) && wait -n
	echo -e "Starting run $job"
	IFS=" " read -r -a temp_array <<< $(sed -n "$job"p temp_job_lines.txt)
	srun -n 1 -c "$SLURM_CPUS_PER_TASK" --exact fastqc -t "$SLURM_CPUS_PER_TASK" -o fastqc "${temp_array[@]}" &
	((++curr_jobs))
done
wait

rm temp_job_lines.txt


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)
echo -e "\nFinished FastQC analysis at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours"
