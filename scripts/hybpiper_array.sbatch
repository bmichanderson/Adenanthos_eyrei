#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Jun 2024
# Description: run an array of HybPiper jobs on Setonix
# Note: requires a Singularity container with HybPiper installed
#	Arguments:
#	-a	additional arguments in " " (optional)
#	-f	target file
#	-r	directory of reads to assemble
#	-s	samples file with one sample per line (reads named: sample_*.fastq.gz)
#	-t	type of target ("dna" default, or "aa")
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=hybpiper
#SBATCH --array=0-156		# this needs to match how many samples to assemble and start from zero
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G		# per job; roughly 1.8 GB per cpu on Setonix
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=hybpiper_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID}"


# set the container
container="/software/projects/pawsey0220/anderson/singularity-containers/hybpiper.sif"


# load module (may need to be updated depending on Setonix updates)
module load singularity/4.1.0-slurm


# parse the command line
while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-a)
	add_args="$2"
	shift
	shift
	;;
	-f)
	target_file="$(readlink -f $2)"
	shift
	shift
	;;
	-r)
	reads_dir="$(readlink -f $2)"
	shift
	shift
	;;
	-s)
	samples_file="$(readlink -f $2)"
	shift
	shift
	;;
	-t)
	target_type="$2"
	shift
	shift
	;;
esac
done

# check args
# create the argument array and append additional arguments if present
arg_array=()
arg_array+=("--diamond")
arg_array+=("--merged")
if [ -n "$add_args" ]; then
	IFS=" " read -r -a temp_array <<< "$add_args"
	for element in "${temp_array[@]}"
	do
		arg_array+=("$element")
	done
	unset IFS
fi

if [ -z "$target_file" ]; then
	echo -e "\nPlease specify a target file\n"
	exit 1
fi

if [ -z "$reads_dir" ]; then
	echo -e "\nPlease specify a reads directory\n"
	exit 1
fi

if [ -z "$samples_file" ]; then
	echo -e "\nPlease specify a samples file\n"
	exit 1
fi

if [ -z "$target_type" ]; then
	targ_piece="-t_dna"
elif [ "$target_type" = "aa" ]; then
	targ_piece="-t_aa"
else
	targ_piece="-t_dna"
fi


# parse the samples file to determine which sample to run
IFS=$'\n' read -d '' -r -a sample_array < $samples_file
sample="${sample_array[$SLURM_ARRAY_TASK_ID]}"
echo "This job will assemble sample $sample"
unset IFS


# run HybPiper

# record the current directory as output directory
outdir="$(pwd)"

# make required directories (if not made yet)
if [ ! -d logs ]; then
	mkdir -p logs
fi
if [ ! -d results ]; then
	mkdir -p results
fi

# make an output directory for this sample
if [ ! -d results/"$sample" ]; then
	mkdir results/"$sample"
fi

# make a working directory and change there
workdir="$(pwd)"/working_"$sample"
if [ ! -d $workdir ]; then
	mkdir -p $workdir
fi
cd $workdir

# run assembly
start="$(date +%s)"

srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
	hybpiper assemble "$targ_piece" "$target_file" -r "$reads_dir"/"$sample"_* \
	--prefix "$sample" --cpu "$SLURM_CPUS_PER_TASK" "${arg_array[@]}" \
	>> "$outdir"/logs/"$sample".out 2>> "$outdir"/logs/"$sample".err

echo "$sample" > "$sample"_done.txt

# collect stats
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
	hybpiper stats "$targ_piece" "$target_file" gene "$sample"_done.txt \
	--seq_lengths_filename "$outdir"/results/"$sample"/"$sample"_lengths \
	--stats_filename "$outdir"/results/"$sample"/"$sample"_stats \
	>> "$outdir"/logs/"$sample".out 2>> "$outdir"/logs/"$sample".err

# collect sequences
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
	hybpiper retrieve_sequences "$targ_piece" "$target_file" dna \
	--sample_names "$sample"_done.txt --fasta_dir "$outdir"/results/"$sample"/dna_seqs \
	>> "$outdir"/logs/"$sample".out 2>> "$outdir"/logs/"$sample".err

srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
	hybpiper retrieve_sequences "$targ_piece" "$target_file" supercontig \
	--sample_names "$sample"_done.txt --fasta_dir "$outdir"/results/"$sample"/supercontig_seqs \
	>> "$outdir"/logs/"$sample".out 2>> "$outdir"/logs/"$sample".err

# collect paralogs
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
	hybpiper paralog_retriever "$targ_piece" "$target_file" "$sample"_done.txt \
	--fasta_dir_all	"$outdir"/results/"$sample"/paralogs_all \
	--fasta_dir_no_chimeras "$outdir"/results/"$sample"/paralogs_no_chimeras \
	--paralog_report_filename "$outdir"/results/"$sample"/paralog_report \
	>> "$outdir"/logs/"$sample".out 2>> "$outdir"/logs/"$sample".err

# remove extraneous lines from the log file
sed -i '/ETA:/d' "$outdir"/logs/"$sample".err
sed -i '/cannot open \/dev\/tty/d' "$outdir"/logs/"$sample".err

# report completion
echo "$sample" >> "$outdir"/done.txt

end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
echo -e "The assembly for sample $sample finished in $duration_mins minutes"

# remove unnecessary files and compress results
mv "$sample"/genes_with_seqs.txt "$sample"/temp
rm "$sample"/*.*
rm "$sample"/*/*.fastq
rm "$sample"/*/"$sample"/*.*
rm -r "$sample"/*/"$sample"/sequences/FAA
mv "$sample"/temp "$sample"/genes_with_seqs.txt
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" \
	tar -c --use-compress-program=pigz -f "$outdir"/results/"$sample"/"$sample".tar.gz "$sample"
cd ..

# remove the working directory
find -P "$workdir" -type f -print0 -o -type l -print0 | xargs -0 munlink
find -P "$workdir" -type d -empty -delete
