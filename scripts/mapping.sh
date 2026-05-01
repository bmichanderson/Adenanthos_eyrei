#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Modified: Mar 2025 (removed module load)
# Description: run paired-end read mapping to a reference on Setonix
# Note: this script should be placed in the launch directory and called by an array script (mapping_array.sbatch)
#	When called, it requires at least three arguments:
#		-1 forward read file
#		-2 reverse read file
#		-r reference fasta
#		-o output directory (default: current directory)
######################


# define software locations
bbmap_loc="/software/projects/pawsey0220/anderson/software/bbmap"
container="/software/projects/pawsey0220/anderson/singularity-containers/samtools.sif"


# set max memory to use (85% of total for Java)
max_mem=$(echo "${SLURM_MEM_PER_NODE}*85/100" | bc)m
echo "This job is going to use $SLURM_CPUS_PER_TASK cores with max mem of $max_mem"


# load modules (may need to be updated depending on Setonix updates)
module load singularity/4.1.0-slurm


# bbmap args
minid="0.9"		# minimum identity for read mapping
ambig="random"	# the way to treat ambiguous mapping; random = randomly assign read; best = first best location
pairlen="400"	# maximum allowed distance between paired reads; insert = pairlen + read1 + read2
pairedonly="t"	# whether to map paired reads only
maxindel="400"	# how large an indel to search for
strictmaxindel="t"	# don't allow indels larger than maxindel


# parse the command line
if [ $# -eq 0 ]; then		# if there are no command arguments
	exit 1
fi
while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-1)
	read1="$(readlink -f $2)"
	shift
	shift
	;;
	-2)
	read2="$(readlink -f $2)"
	shift
	shift
	;;
	-r)
	ref="$(readlink -f $2)"
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
if [ -z "$read1" ] || [ -z "$read2" ] || [ -z "$ref" ]; then		# if input files are not specified
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
file_name=$(basename ${read1})
prefix=${file_name/_*/}
logfile=mapping_"$prefix".log


# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting read mapping at $(date)\n**********" >> $logfile 2>&1


# Run read mapping
"$bbmap_loc"/bbmap.sh ref="$ref" in1="$read1" in2="$read2" outm=map.sam bs=bs.sh nodisk \
	minid="$minid" ambiguous="$ambig" pairedonly="$pairedonly" pairlen="$pairlen" \
	maxindel="$maxindel" strictmaxindel="$strictmaxindel" \
	threads="$SLURM_CPUS_PER_TASK" -Xmx"$max_mem" >> $logfile 2>&1
echo -e "\n\n********\nFinished mapping\n*********\n" >> $logfile 2>&1
singularity exec -H "$(pwd)" "$container" bash bs.sh >> $logfile 2>&1
rm bs.sh map.sam


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)
echo -e "\nFinished read mapping at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours\n" >> $logfile 2>&1

echo "Job ${SLURM_ARRAY_TASK_ID} for sample "$prefix" finished in $duration seconds or $duration_mins minutes"
