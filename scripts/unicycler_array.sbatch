#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Description: run an array of Unicycler assembly jobs on Setonix, with optional read mapping to a reference
# Note: this requires a file with sample names for mapping and (optionally) reference fasta file(s)
#	to map to (tab separated, one per line)
#	i.e.	sample1 \t /path/to/ref1
#			sample2 \t /path/to/ref2
#			...
#	or		sample1
#			sample2
#			...
#	It expects reads in the current directory named in the format: `sample_...R1.fastq.gz` and `sample_...R2.fastq.gz`
#	It needs a `unicycler.sh` executable file in the current directory
#	It will create a folder for each sample to be mapped (if they don't exist yet)
# Note: set the "--array=__-__" variable to the number of samples and starting from 0!
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=unicycler
#SBATCH --array=0-9		# this needs to match how many samples will be assembled
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=unicycler_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID}"


# check arguments
if [ -z $1 ]; then
	echo -e "Please provide a text file with names of samples and (optionally) reference fastas\n"
	exit 1
fi


# Determine the sample and (optionally) reference
# first count columns
ncols=$(head -n 1 $1 | awk '{print NF}')
# then assign vars
samples=($(cut -f 1 $1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}
if [ $ncols -gt 1 ]; then
	refs=($(cut -f 2 $1))
	ref=${refs[$SLURM_ARRAY_TASK_ID]}
fi


# Run the Unicycler script for the given pair of read files and (optionally) reference
echo "This job is for sample ${sample}"
if [ $ncols -gt 1 ]; then
	srun -N 1 -n 1 -c ${SLURM_CPUS_PER_TASK} ./unicycler.sh -1 ${sample}_*R1.fastq.gz \
		-2 ${sample}_*R2.fastq.gz -r ${ref} -o ${sample}
else
	srun -N 1 -n 1 -c ${SLURM_CPUS_PER_TASK} ./unicycler.sh -1 ${sample}_*R1.fastq.gz \
		-2 ${sample}_*R2.fastq.gz -o ${sample}
fi
