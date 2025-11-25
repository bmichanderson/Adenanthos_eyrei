#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Description: run an array of read mapping jobs on Setonix
# Note: this requires a file with sample names for mapping and reference fasta file(s) to map to (tab separated, one per line)
#	i.e.	sample1 \t /path/to/ref1
#			sample2 \t /path/to/ref2
#			...
#	It expects reads in the current directory named in the format: `sample_...R1.fastq.gz` and `sample_...R2.fastq.gz`
#	It needs a `mapping.sh` executable file in the current directory
#	It will create a folder for each sample to be mapped (if they don't exist yet)
# Note: set the "--array=__-__" variable to the number of samples and starting from 0!
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=mapping
#SBATCH --array=0-125		# this needs to match how many samples will be mapped
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --output=map_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID}"


# check arguments
if [ -z $1 ]; then
	echo -e "Please provide a text file with names of samples and reference fastas\n"
	exit 1
fi


# Determine the sample and reference (based on first argument text file)
samples=($(cut -f 1 $1))
refs=($(cut -f 2 $1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}
ref=${refs[$SLURM_ARRAY_TASK_ID]}


# Run the mapping script for the given pair of read files and reference
echo "This job is for sample ${sample}"
srun -N 1 -n 1 -c ${SLURM_CPUS_PER_TASK} ./mapping.sh -1 ${sample}_*R1.fastq.gz \
	-2 ${sample}_*R2.fastq.gz -r ${ref} -o ${sample}
