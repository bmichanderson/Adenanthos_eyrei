#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Modified: Mar 2025 (cleaned up)
# Description: run an array of GetOrganelle assemblies from SPAdes assemblies on Setonix
# Note: a Singularity container with SPAdes installed ("getorganelle.sif") is required
#	An executable "plastome_from_assembly.sh" script is required in the current directory
#	(see that script for required databases to download)
#	Each job will need input (SPAdes assembly) and output directories specified
#	To do this, the script takes a text file (arg1) with folder names/paths, one per line
#	It assumes that the folder (base)names correspond to samples
#	It will make a directory {sample} for the assembly in the launch directory, or use
#	an existing folder if that folder is already there
# Note: set the "--array=__-__" variable to the number of samples and starting from 0!
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=getorganelle
#SBATCH --array=0-125		# this needs to match how many samples to assemble
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G		# roughly 1.8 GB per core on Setonix
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=getorg_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID}"


# check argument
if [ -z $1 ]; then
	echo -e "Please provide a text file with paths to SPAdes assembly folders\n"
	exit 1
fi


# Determine the folder paths (based on first argument text file)
folders=($(cat $1))
folder=${folders[$SLURM_ARRAY_TASK_ID]}
folder=$(readlink -f $folder)


# Run the spades script for the given sample
srun -N 1 -n 1 -c ${SLURM_CPUS_PER_TASK} ./plastome_from_assembly.sh -s $folder \
	-o $(basename $folder)
