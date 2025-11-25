#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Modified: Mar 2025 (added -e argument to allow for exact bp extension specification)
# Description: run an array of Samtools summary jobs on Setonix
# Note: The main agrument (-f) is a file with sample names, which should match folders in
#	the launch directory that contain a single BAM mapping file (ideally, created with `mapping.sh`)  
#	It needs a `samtools_summary.sh` executable file in the current directory
#	If wanting to run on an extended reference, specify the extension (-e) in bp
# Note: set the "--array=__-__" variable to the number of samples and starting from 0!
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=samtools
#SBATCH --array=0-75		# this needs to match how many samples will be run
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --output=samtools_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID}"


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
	-f)
	folders_file="$(readlink -f $2)"
	shift
	shift
	;;
esac
done


# check arguments
if [ -z "$folders_file" ]; then
	echo -e "Please provide a text file with names of samples matching folders in the launch directory\n"
	exit 1
fi
if [ -z "$extend" ]; then
	extend=0
fi


# Determine the sample 
samples=($(cut -f 1 $folders_file))
sample=${samples[$SLURM_ARRAY_TASK_ID]}


# Run the summary script for the given sample
echo "This job is for sample ${sample}"
srun -N 1 -n 1 -c ${SLURM_CPUS_PER_TASK} ./samtools_summary.sh -s $sample -e $extend
