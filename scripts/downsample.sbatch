#!/bin/bash -l

######################
# Author: B.M. Anderson
# Date: Oct 2023
# Modified: Mar 2025 (file naming; added options for approximate number output and read length priority)
# Description: downsample Illumina paired-end reads on Setonix
# Note:	specify the read files as args -1 and -2;
#	if wanting to downsample to a proportion, use -p {proportion} [default 0.25], or
#	if wanting to downsample to an approximate number of reads, use -a {number};
#	if wanting to prioritize longer reads, use -r "yes" (can only be used when setting an exact number);
#	The resulting files will be output in the launch directory
# Note: for file output naming, this script assumes reads are named with a prefix before an underscore
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=downsample
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --output=downsample_%j.log


# set max memory to use (85% of total for Java)
max_mem=$(echo "${SLURM_MEM_PER_NODE}*85/100" | bc)m
echo "This job is going to use $SLURM_CPUS_PER_TASK cores with max mem of $max_mem"


# define software location
bbmap_loc="/software/projects/pawsey0220/anderson/software/bbmap"


# set call args
call_args=()


# parse the command line
approx_present="False"
proportion_present="False"
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
	-a)
	approx="$2"
	approx_present="True"
	shift
	shift
	;;
	-p)
	proportion="$2"
	proportion_present="True"
	shift
	shift
	;;
	-r)
	read_prio="$2"
	shift
	shift
	;;
esac
done

# check command args
if [ -z "$read1" ] || [ -z "$read2" ]; then
	echo -e "\nMissing read files!\n"
	exit 1
fi
if [ "$approx_present" == "True" -a "$proportion_present" == "True" ]; then
	echo -e "\nSpecify either exact (-e) or proportion (-p), but not both\n"
	exit 1
elif [ "$approx_present" == "True" ]; then
	call_args+=("samplereadstarget=${approx}")
	if [ "$read_prio" == "yes" ]; then
		call_args+=("prioritizelength=t")
	fi
elif [ "$proportion_present" == "True" ]; then
	call_args+=("samplerate=${proportion}")
else
	call_args+=("samplerate=0.25")
fi


# set the output prefix (assuming format: sampleID_*.fastq.gz)
file_name=$(basename $read1)
prefix=${file_name/_*/}


# downsample
echo -e "\n*****\nDownsampling\n*****\n"
"$bbmap_loc"/reformat.sh in1="$read1" in2="$read2" out1="$prefix"_down_R#.fastq.gz \
	${call_args[@]} -Xmx"$max_mem"
echo -e "\n\n"


# run readlength stats on the final files
"$bbmap_loc"/readlength.sh in1="$prefix"_down_R1.fastq.gz in2="$prefix"_down_R2.fastq.gz
