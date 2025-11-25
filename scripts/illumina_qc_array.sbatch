#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Jun 2024
# Description: run an array of Illumina QC jobs on Setonix
# Note: this will run on all *.gz files in the directory where it is launched
#	It expects these to be named in the format: `prefix_...R1.fastq.gz` and `prefix_...R2.fastq.gz`
#	It needs an `illumina_qc.sh` executable file in the current directory
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=illuminaqc
#SBATCH --array=0-125		# this needs to match how many pairs of read files are in the current directory
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32		# this can be increased even more if desired, since memory is already requesting a full node (128 cpus)
#SBATCH --mem=230G			# 230G is the max for a work node
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --output=illuminaqc_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID}"


# Determine the file prefixes (based on all the gzipped files in the current directory)
fileprefixes=($(ls *.gz | cut -f 1 -d "_" | uniq))
fileprefix=${fileprefixes[$SLURM_ARRAY_TASK_ID]}


# Run the QC script for the given pair of read files
echo "This job is for sample ${fileprefix}"
srun -N 1 -n 1 -c ${SLURM_CPUS_PER_TASK} ./illumina_qc.sh -1 ${fileprefix}*R1.fastq.gz \
	-2 ${fileprefix}*R2.fastq.gz
