#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Modified: Mar 2025 (adjusted reporting and removed one module)
# Description: run a GetOrganelle assembly of a plastome from an input SPAdes assembly on Setonix
# Note: a Singularity container with the dependencies ("getorganelle.sif") is required
#	The script requires a path to a folder with a SPAdes assembly (includes "scaffolds.paths" and "assembly_graph.fastg")
#	The output will be in a folder "getorganelle_out" in the working directory, but the working directory
#	can be specified with the second argument
#	This script is meant to be called by an array script
#		plastome_from_assembly.sh -s <spades_folder> -o <outdir>
#
#	Pre-downloaded databases are needed in the launch directory (/.GetOrganelle)
#	Download them interactively like so:
#	salloc -p copy -n 1 -N 1 -c 4 --mem=8G -A pawsey0220 --time=01:00:00
#		module load singularity/4.1.0-slurm
#		container="/software/projects/pawsey0220/anderson/singularity-containers/getorganelle.sif"
#		singularity exec -H "$(pwd)" "$container" get_organelle_config.py \
#			--config-dir "$(pwd)"/.GetOrganelle --add embplant_pt,embplant_mt
######################


# define container to use
container="/software/projects/pawsey0220/anderson/singularity-containers/getorganelle.sif"


# load module (may need to be updated depending on Setonix updates)
module load singularity/4.1.0-slurm


# Set launch directory variable (where the database folder should be)
launchdir="$(pwd)"


# parse the command line
if [ $# -eq 0 ]; then		# if there are no command arguments
	exit 1
fi
while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-s)
	folder="$(readlink -f $2)"
	shift
	shift
	;;
	-o)
	out_dir="$(readlink -f $2)"
	shift
	shift
	;;
esac
done

# check args
if [ -z "$folder" ]; then		# if input folder not specified
	exit 1
fi
if [ -z "$out_dir" ]; then		# if output path is not specified
	out_dir="$(pwd)"
fi


# change to out_dir (creating it if it doesn't exist)
if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi
cd "$out_dir"


# report how many cores are going to be used and what sample
echo "Job ${SLURM_ARRAY_TASK_ID} will assemble sample $(basename $folder) using $SLURM_CPUS_PER_TASK cores"


# Set up a log file
logfile=getorganelle_"$SLURM_JOB_ID".log


# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting plastome assembly at $(date)\n**********\n" >> $logfile 2>&1


# Run the assembly
singularity exec -H "$(pwd)" "$container" get_organelle_from_assembly.py \
	--config-dir "$launchdir"/.GetOrganelle -g "$folder"/assembly_graph.fastg \
	-o getorganelle_out --spades-out-dir "$folder" -F embplant_pt \
	-t "$SLURM_CPUS_PER_TASK" >> $logfile 2>&1


# Remove unnecessary files
rm getorganelle_out/initial_assembly_graph.fastg


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)
echo -e "\nFinished plastome assembly at $(date) after running for" \
	"$duration_mins minutes or $duration_hours hours" >> $logfile 2>&1
echo "Job ${SLURM_ARRAY_TASK_ID} for sample $(basename $folder) finished in $duration_mins minutes"
