#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: June 2024
# Modified: May 2025 (slightly adjusted execution, reporting and settings)
# Description: run alignment, clean up and phylogenetic analyses for a directory of fasta loci
# Note:	requires a Singularity container with all the dependencies/programs (phylo.sif)
#	Will output results in the launch directory and will create temp folders and a folder called "in_align"
#	Arguments:
#	-a analyses to run: "all" [default] "concat" "loci" "none"
#	-c clean alignments: "y" [default] or "n" (default removing sites with < 50% data and samples with < 75% data)
#	-con concordance analysis: "y" [default] or "n"
#	-f fasta directory
#	-l likelihood mapping: "y" or "n" [default]
#	-p polytomy simplification: "y" or "n" [default] (collapse branches in locus trees with < 50 UFB)
#	-r re-alignment: "y" [default] or "n" (re-align the input fasta files with MAFFT)
#
#	If wanting to run tasks in parallel (e.g. both concat and loci, or likelihood mapping), launch with multiple tasks, e.g.
#	sbatch --ntasks=5 --cpus-per-task=12 align_phylo.sbatch -a "all" -l "y"
#	sbatch --ntasks=2 --cpus-per-task=32 align_phylo.sbatch -a "all"
#
#	For better efficiency for many loci, run analyses separately, e.g.
#	sbatch --ntasks=1 --cpus-per-task=8 align_phylo.sbatch -a "concat"
#	sbatch --ntasks=1 --cpus-per-task=64 align_phylo.sbatch -a "loci"
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=align_phylo
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --export=NONE
#SBATCH --output=align_phylo_%j.log


# Set the container and cleaning script
container="/software/projects/pawsey0220/anderson/singularity-containers/phylo.sif"
clean_script="/software/projects/pawsey0220/anderson/scripts/clean_alignment.py"


# Load modules
module load python/3.11.6 singularity/4.1.0-slurm


# Settings
max_mem="32G"	# ASTRAL Java memory setting
polytomy="50"	# percentage threshold for collapsing low support branches in locus tress
max_miss="50"	# maximum percentage missing or N per site to keep (clean step)
sample_miss="75"	# maximum percentage missing or N per sample to keep (clean step)


# Set up a master log file
logfile=master_"${SLURM_JOBID}".log
logfile="$(readlink -f $logfile)"


# Parse the command line
while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-a)
	analysis="$2"
	shift
	shift
	;;
	-c)
	clean="$2"
	shift
	shift
	;;
	-con)
	concord="$2"
	shift
	shift
	;;
	-f)
	align_dir="$(readlink -f $2)"
	shift
	shift
	;;
	-l)
	likemap="$2"
	shift
	shift
	;;
	-p)
	poly="$2"
	shift
	shift
	;;
	-r)
	realign="$2"
	shift
	shift
	;;
esac
done

# Check args
if [ -z "$analysis" ]; then
	analysis="all"
elif [ "$analysis" == "all" ] || [ "$analysis" == "concat" ] || [ "$analysis" == "loci" ] || [ "$analysis" == "none" ]; then
	analysis="$analysis"
else
	echo -e "\nPlease specify analysis as \"all\", \"concat\", \"loci\" or \"none\"\n"
	exit 1
fi

if [ -z "$clean" ]; then
	clean="y"
fi

if [ -z "$concord" ]; then
	concord="y"
fi

if [ -z "$align_dir" ]; then
	echo -e "\nPlease specify where the fasta alignments are\n"
	exit 1
fi

if [ -z "$likemap" ]; then
	likemap="n"
fi

if [ -z "$poly" ]; then
	poly="n"
fi

if [ -z "$realign" ]; then
	realign="y"
fi


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting alignment clean and phylogenetics run at $(date)\n**********" | tee -a $logfile 2>&1


# Make a directory for the alignments to be put in and to be input to phylogenetic analysis
in_dir="in_align"
if [ ! -d "$in_dir" ]; then
	mkdir -p "$in_dir"
fi
in_dir=$(readlink -f $in_dir)


# Re-align if not turned off
if [ "$realign" == "y" ]; then
	echo -e "\n****\tRunning initial re-alignments" >> $logfile 2>&1
	align_start="$(date +%s)"
	# make a temp directory and re-align
	if [ ! -d temp_align ]; then
		mkdir temp_align
	fi

	for align in "$align_dir"/*.fasta
	do
		align_name="$(basename $align)"
		srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
			mafft --quiet --auto --thread "$SLURM_CPUS_PER_TASK" "$align" > temp_align/"$align_name"
	done
	align_end="$(date +%s)"
	duration="$(( $align_end - $align_start ))"
	echo -e "\tFinished re-alignments in $duration seconds\n" >> $logfile 2>&1
	align_dir=temp_align
else
	echo -e "\n****\tSkipping re-alignment\n" >> $logfile 2>&1
fi


# Clean the alignments if not turned off
if [ "$clean" == "y" ]; then
	echo -e "\n****\tCleaning alignments" >> $logfile 2>&1
	clean_start="$(date +%s)"
	if [ ! -d temp_clean ]; then
		mkdir temp_clean
	fi
	cd temp_clean
	ln -s ../"$align_dir"/*.fasta .
	for align in *.fasta
	do
		echo -e "Cleaning alignment $align"
		python3 "$clean_script" -p "$max_miss" -s "$sample_miss" "$align"
	done
	mv *_clean.fasta "$in_dir"
	cd ..
	clean_end="$(date +%s)"
	duration="$(( $clean_end - $clean_start ))"
	echo -e "\tFinished cleaning alignments in $duration seconds\n" >> $logfile 2>&1
else
	echo -e "\n****\tSkipping cleaning step\n" >> $logfile 2>&1
	rsync -rutL "$align_dir"/*.fasta "$in_dir"
fi


# Remove the temporary alignment and cleaning folders, if present
if [ -d temp_align ]; then
	rm -r temp_align
fi
if [ -d temp_clean ]; then
	rm -r temp_clean
fi


# Run a concatenation tree using IQTREE (in the background if also running loci and multiple tasks available)
if [ "$analysis" == "all" ]; then
	# if the job has been launched with multiple tasks, run concatenation in the background
	if [ "$SLURM_NTASKS" -gt 1 ]; then
		echo -e "\n****\tRunning ML analysis of the concatenated alignment in the background (using maximum $SLURM_CPUS_PER_TASK CPUs)\n" >> $logfile 2>&1
		srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
			iqtree --threads-max "$SLURM_CPUS_PER_TASK" -T AUTO --ufboot 1000 --sampling GENESITE \
			-m MFP --prefix concat -p "$in_dir" > /dev/null 2>&1  &
	else
		# run in serial
		echo -e "\n****\tRunning ML anlaysis of the concatenated alignment (using maximum $SLURM_CPUS_PER_TASK CPUs)\n" >> $logfile 2>&1
		srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
			iqtree --threads-max "$SLURM_CPUS_PER_TASK" -T AUTO --ufboot 1000 --sampling GENESITE \
			-m MFP --prefix concat -p "$in_dir" > /dev/null 2>&1
	fi
elif [ "$analysis" == "concat" ]; then
	# run the analysis
	echo -e "\n****\tRunning ML anlaysis of the concatenated alignment (using maximum $SLURM_CPUS_PER_TASK CPUs)\n" >> $logfile 2>&1
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
		iqtree --threads-max "$SLURM_CPUS_PER_TASK" -T AUTO --ufboot 1000 --sampling GENESITE \
		-m MFP --prefix concat -p "$in_dir" > /dev/null 2>&1
else
	echo -e "\n****\tNot running analysis of the concatenated alignment\n" >> $logfile 2>&1
fi


# Run trees for each alignment
if [ "$analysis" == "all" ] || [ "$analysis" == "loci" ]; then
	echo -e "\n****\tRunning ML analyses of loci alignments (using maximum $SLURM_CPUS_PER_TASK CPUs)\n" >> $logfile 2>&1
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
		iqtree --threads-max "$SLURM_CPUS_PER_TASK" -T AUTO --ufboot 1000 -m MFP \
		--prefix loci -S "$in_dir" > /dev/null 2>&1
else
	echo -e "\n****\tNot running analyses of individual loci alignments\n" >> $logfile 2>&1
fi


# Run likelihood mapping if requested
if [ "$likemap" == "y" ]; then
	mkdir likemap
	ln -s "$in_dir"/*.fasta likemap/
	cd likemap
	echo -e "\n****\tRunning likelihood mapping for loci alignments (using maximum $SLURM_CPUS_PER_TASK CPUs per locus)\n" >> $logfile 2>&1
	if [ "$analysis" == "all" ] && [ "$SLURM_NTASKS" -gt 1 ]; then
		max_jobs=$((SLURM_NTASKS - 1))
	else
		max_jobs="$SLURM_NTASKS"
	fi
	curr_jobs=0
	for alignment in *.fasta
	do
		((curr_jobs >= max_jobs)) && wait -n
		srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" singularity exec -H "$(pwd)" "$container" \
			iqtree -T "$SLURM_CPUS_PER_TASK" -s "$alignment" --lmap 10000 --quartetlh -n 0 \
			--prefix lm_"${alignment/.fasta/}" > /dev/null 2>&1 &
		((++curr_jobs))
	done
	wait
	rm *.fasta *.phy *.log *.gz *.treefile *.eps *.svg
	cd ..
else
	echo -e "\n****\tNot running likelihood mapping\n" >> $logfile 2>&1
fi


# Wait for concatenation to finish if running in the background
if [ "$analysis" == "all" ] && [ "$SLURM_NTASKS" -gt 1 ]; then
	wait
	echo -e "\tFinished background job for the ML analysis of the concatenated alignment" >> $logfile 2>&1
fi


# Reduce low-support branches to polytomies, if requested
if [ "$poly" == "y" ]; then
	if [ -f "loci.treefile" ]; then
		echo -e "\n****\tCollapsing branches with < $polytomy UFB to polytomies" >> $logfile 2>&1
		srun -N 1 -n 1 -c 1 --exact singularity exec -H "$(pwd)" "$container" \
			nw_ed loci.treefile "i & b < $polytomy" o > loci_poly.tre
		mv loci_poly.tre loci.treefile
	else
		echo -e "\n****\tNo loci tree file found, so not collapsing branches" >> $logfile 2>&1
	fi
fi


# Calculate concordance factors on the concatenation tree if requested and if the relevant files exist
if [ "$concord" == "y" ]; then
	if [ -f "concat.treefile" ]; then
		# to re-use the concatenation model and tree, need to copy the checkpoint and
		# model files to match the new prefix
		cp concat.model.gz concord_scf.model.gz
		cp concat.ckp.gz concord_scf.ckp.gz
		if [ -f "loci.treefile" ]; then
			echo -e "\n****\tCalculating gene and site concordance factors" >> $logfile 2>&1
			srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
				iqtree -T "$SLURM_CPUS_PER_TASK" -t concat.treefile --gcf loci.treefile \
				--prefix concord_gcf > /dev/null 2>&1
			srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
				iqtree -T "$SLURM_CPUS_PER_TASK" -te concat.treefile --scfl 10000 \
				--prefix concord_scf -p "$in_dir" -blfix --undo > /dev/null 2>&1
		else
			echo -e "\n****\tCalculating site concordance factors" >> $logfile 2>&1
			srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
				iqtree -T "$SLURM_CPUS_PER_TASK" -te concat.treefile --scfl 10000 \
				--prefix concord_scf -p "$in_dir" -blfix --undo > /dev/null 2>&1
		fi
	fi
else
	echo -e "\n****\tNot running concordance factor analyses\n" >> $logfile 2>&1
fi


# Run an ASTRAL tree using the gene trees (if they exist)
if [ -f "loci.treefile" ]; then
	echo -e "\n****\tRunning ASTRAL analyses on the loci treefile" >> $logfile 2>&1
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
		java -Xmx"$max_mem" -jar /Astral/astral.5.7.1.jar -i loci.treefile \
		-o astral.tre --branch-annotate 2 > astral_screenlog.txt 2>&1
	# add a polytomy test that will put p-values on branches
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$container" \
		java -Xmx"$max_mem" -jar /Astral/astral.5.7.1.jar -i loci.treefile \
		-o astral_poly.tre --branch-annotate 10 >> astral_screenlog.txt 2>&1
else
	echo -e "\n****\tNo loci tree file found, so not running ASTRAL" >> $logfile 2>&1
fi


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished alignment clean and phylogenetics run at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" | tee -a $logfile 2>&1
