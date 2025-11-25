#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Description: run an array of SPAdes assemblies on Setonix
# Note: a Singularity container with SPAdes installed ("getorganelle.sif") is required
#	An executable "spades.sh" script is required in the current directory
#	Each job will need two input read files and an output directory specified
#	To do this, the script assumes all the *.gz files in the current directory are input
#	and each is named in the pattern "sample_*R1.fastq.gz" and "sample_*R2.fastq.gz"
#	It will make a directory "sample" for the assembly
# Note: set the "--array=__-__" variable to the number of samples and starting from 0!
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=spades
#SBATCH --array=0-125		# this needs to match how many samples to assemble
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G		# roughly 1.8 GB per cpu on Setonix
#SBATCH --time=12:00:00
#SBATCH --export=NONE
#SBATCH --output=spades_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID}"


# Determine the file prefixes (based on all the gzipped files in the current directory)
fileprefixes=($(ls *.gz | cut -f 1 -d "_" | uniq))
fileprefix=${fileprefixes[$SLURM_ARRAY_TASK_ID]}


# Run the spades script for the given sample
echo "This job is for sample ${fileprefix}"
srun -N 1 -n 1 -c ${SLURM_CPUS_PER_TASK} --mem ${SLURM_MEM_PER_NODE} ./spades.sh \
	-1 ${fileprefix}_*R1.fastq.gz -2 ${fileprefix}_*R2.fastq.gz -o ${fileprefix}
