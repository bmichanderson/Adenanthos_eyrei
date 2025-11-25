#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Modified: Mar 2025 (updated modules loaded)
# Description: run a single SPAdes assembly on Setonix
# Note: a Singularity container with SPAdes installed ("getorganelle.sif") is required
#	This script should be placed in the working directory and called by the array script ("spades_array.sbatch")
######################


# Set OMP variable for SPAdes multithreading
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"


# define container to use
container="/software/projects/pawsey0220/anderson/singularity-containers/getorganelle.sif"


# set max memory to use (SPAdes wants an integer in GB; set to 90% to be conservative)
max_mem=$(echo "${SLURM_MEM_PER_NODE}*90/100000" | bc)
echo "This job is going to use $SLURM_CPUS_PER_TASK cores with max mem of ${max_mem}GB"


# load module (may need to be updated depending on Setonix updates)
module load singularity/4.1.0-slurm


# parse the command line
if [ $# -eq 0 ]; then		# if there are no command arguments
	exit 1
fi

while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-1)
	reads1="$(readlink -f $2)"
	shift
	shift
	;;
	-2)
	reads2="$(readlink -f $2)"
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
if [ -z "$reads1" ] || [ -z "$reads2" ]; then		# if input files are not specified
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


# make a log file based on the read name
file_name=$(basename ${reads1})
prefix=${file_name/_*/}
logfile=spades_"$prefix".log


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting SPAdes assembly at $(date)\n**********" >> $logfile 2>&1


# Run assembly
singularity exec -H "$(pwd)" "$container" spades.py -1 "$reads1" -2 "$reads2" \
	-o spades_out --only-assembler -t "$SLURM_CPUS_PER_TASK" --memory "$max_mem" \
	-k 21,45,65,85,105,125 >> $logfile 2>&1


# Remove extra spades output, keeping the scaffolds and assembly graphs
mv "$(pwd)"/spades_out/scaffolds.fasta .
mv "$(pwd)"/spades_out/scaffolds.paths .
mv "$(pwd)"/spades_out/assembly_graph_with_scaffolds.gfa .
mv "$(pwd)"/spades_out/assembly_graph.fastg .
rm -r spades_out/


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)
echo -e "\nFinished SPAdes assembly at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours\n" >> $logfile 2>&1
echo "Job ${SLURM_ARRAY_TASK_ID} for sample "$prefix" finished in $duration_mins minutes"
