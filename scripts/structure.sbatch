#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Sep 2021
# Modified: Oct 2021 (allow specifying starting rep number); Nov 2021 (allow specifying infile; cleaned up);
#	Aug 2022 (move to Setonix; increase the number of tasks to run in parallel; now use a container); Mar 2025 (updated module version)
#	Feb 2026 (output log files as run progresses)
# Description: run a Structure analysis
# Note:	ensure ntasks is set to the number of reps
# Note:	call the script with a specified K value (-k), working directory (-w) [default current directory]
#	with files "mainparams" and "extraparams", and the Structure input file (-s), which can be omitted
#	if "infile" exists in the working directory
# Note:	if breaking the job up more, the starting rep can be specified (-r), but this shouldn't be needed now with 40 cores  
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=struct
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --export=NONE
#SBATCH --output=struct_%j.log


# set the location of the container
container="/software/projects/${PAWSEY_PROJECT}/anderson/singularity-containers/structure.sif"


# load modules (may need to be updated depending on Setonix updates)
module load singularity/4.1.0-slurm


# Parse arguments
if [ $# -eq 0 ]; then		# if there are no command arguments
	exit 1
fi
while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-k)
	k_val="$2"
	shift
	shift
	;;
	-w)
	workdir="$(readlink -f $2)"
	shift
	shift
	;;
	-r)
	start_rep="$2"
	shift
	shift
	;;
	-s)
	str_file="$(readlink -f $2)"
	shift
	shift
	;;
esac
done

# check input
if [ -z "$k_val" ]; then
	echo -e "Please specify a K value\n"
	exit 1
fi
if [ -z "$workdir" ]; then
	workdir="$(pwd)"
fi
if [ -z "$start_rep" ]; then
	start_rep="1"
fi
if [ -z "$str_file" ]; then
	str_file="infile"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\nStarting Structure analysis at $(date)\n"
echo -e "Working directory is $workdir\n"


# change to the working directory and start the analysis
cd "$workdir"
echo -e "\nRunning a Structure analysis\n"


# run structure for 40 reps
echo -e "\nRunning for K $k_val\n"

for ((rep=${start_rep}; rep<=((${start_rep} + 39)); rep++))
do
	echo -e "\nStarting rep $rep\n"
	srun -N 1 -n 1 --exact singularity run "$container" -m mainparams -e extraparams -K "$k_val" -i "$str_file" \
		-o output_"$k_val"_"$rep" -D "$rep""$k_val""$RANDOM" > log_"$k_val"_"$rep".txt &
done
wait


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)
echo -e "\nFinished Structure analysis at $(date) after running for" \
	"$duration_mins minutes, or $duration_hours hours"
