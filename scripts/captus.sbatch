#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: May 2025
# Description: run the Captus pipeline (https://github.com/edgardomortiz/Captus) on Setonix
# Note: requires a Singularity container with Captus installed
#	Adjust cpus and memory as appropriate
#	Default settings for each step can be adjusted in this script (make a copy)
#	Launch the job from a directory where there is existing output or where output should go
#	Arguments:
#	-s	Step(s) to run: 1 - cleaning, 2 - assembly, 3 - extraction, 4 - alignment
#		If wanting to run multiple steps, specify without spaces, e.g. "-s 23" for steps 2 and 3
#		It doesn't make sense to run later steps without existing folders/files from earlier
#	-r	Read files (raw or cleaned, depending on the starting step)
#		Can be specified as a directory (e.g. "01_clean_reads") or list of files
#	-t	Nuclear targets file (can be of protein, preferably, or DNA sequences)
#		Entries can be named similarly to: ">taxon-name-locus [anything]" with final "-" splitting
#	-m	Mode to run for extraction (step 3): n - nuclear, p - plastid, m - mitochondrial, e - extras
#		Default is to run with only nuclear ("-m n")
#		If wanting multiple modes, specify without spaces, e.g. "-m np" for nuclear and plastid
#		Multiple modes with separate target files is not implemented here (but possible)
#		If multiple selected, will use the Captus default references for plastid/mitochondrial
#		The "e" extras option will attempt to cluster contigs and find new markers
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=captus
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=230G
#SBATCH --time=12:00:00
#SBATCH --export=NONE
#SBATCH --out=captus-%j.log


# set the container
container="/software/projects/pawsey0220/anderson/singularity-containers/captus.sif"


# load module (may need to be updated depending on Setonix updates)
module load singularity/4.1.0-slurm


# set default running parameters for parallelisation
ram=$(echo "${SLURM_MEM_PER_NODE}/1000" | bc)			# SLURM_MEM_PER_NODE is in MB, ram is specified in GB
threads="${SLURM_CPUS_PER_TASK}"
#concurrent=$(echo "${SLURM_CPUS_PER_TASK}/8" | bc)		# change this as desired for how many cpus per sample (8)
concurrent="auto"		# to allow Captus to determine optimal execution


# parse the command line
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
	-s)
	steps_string="$2"
	shift
	shift
	;;
	-r)
	reads="$(readlink -f $2)"
	shift
	shift
	;;
	-t)
	targets_file="$(readlink -f $2)"
	shift
	shift
	;;
	-m)
	mode_string="$2"
	shift
	shift
	;;
esac
done

# check args
run_step1=false
run_step2=false
run_step3=false
run_step4=false
mode_n=false
mode_p=false
mode_m=false
mode_e=false

if [ -z "$steps_string" ]; then		# no steps specified
	echo -e "\nPlease specify which steps to run with \"-s\"\n"
	exit 1
else
	if [[ "$steps_string" == *"1"* ]]; then
		run_step1=true
	fi
	if [[ "$steps_string" == *"2"* ]]; then
		run_step2=true
	fi
	if [[ "$steps_string" == *"3"* ]]; then
		run_step3=true
	fi
	if [[ "$steps_string" == *"4"* ]]; then
		run_step4=true
	fi
fi

if [ -z "$mode_string" ]; then		# no mode specified
	mode_n=true
else
	if [[ "$mode_string" == *"n"* ]]; then
		mode_n=true
	fi
	if [[ "$mode_string" == *"p"* ]]; then
		mode_p=true
	fi
	if [[ "$mode_string" == *"m"* ]]; then
		mode_m=true
	fi
	if [[ "$mode_string" == *"e"* ]]; then
		mode_e=true
	fi
fi


# record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting Captus run at $(date)\n**********"


# run Captus

## step 1 (cleaning)
if [ "$run_step1" = true ]; then
	# defaults
	output="./01_clean_reads/"
	adaptors="ALL"		# other options: "Illumina" "BGI"
	trimq="13"			# PHRED quality to trim ends
	maq="16"			# minimum average quality across read
	ftl="0"				# force trim left

	# run the cleaning
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity run -H "$(pwd)" "$container" \
		clean --reads "$reads" --out "$output" --adaptor_set "$adaptors" --trimq "$trimq" \
		--maq "$maq" --ftl "$ftl" --ram "$ram" --threads "$threads" --concurrent "$concurrent"
fi


## step 2 (assembly)
if [ "$run_step2" = true ]; then
	# defaults
	if [ "$run_step1" = true ]; then		# cleaning was run
		clean_reads="01_clean_reads"
	else
		clean_reads="$reads"
	fi
	output="./02_assemblies/"
	preset="CAPSKIM"		# other options: "RNA" "WGS"
	# note: can also specify kmers, min depth and prune level more exactly, but preset should be fine
	tmpdir="$(pwd)/tmp"		# temp directory for MEGAHIT
	mindepth="1.5"			# min depth for retained contigs

	# make the temp directory if it doesn't exist
	if [ ! -d "$tmpdir" ]; then
		mkdir -p "$tmpdir"
	fi

	# run the assembly
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity run -H "$(pwd)" "$container" \
		assemble --reads "$clean_reads" --out "$output" --preset "$preset" --tmp_dir "$tmpdir" \
		--min_contig_depth "$mindepth" --ram "$ram" --threads "$threads" --concurrent "$concurrent"
fi


## step 3 (extraction)
if [ "$run_step3" = true ]; then
	# defaults
	in_assemblies="02_assemblies"		# if starting here, ensure this folder is present
	output="./03_extractions/"
	if [ -z "$targets_file" ]; then		# if no targets file provided
		targets="Mega353"		# other options: "Angiosperms353"
	else
		targets="$targets_file"
	fi

	# set the ref args depending on mode
	refarray=()
	if [ "$mode_n" = true ]; then
		refarray+=("--nuc_refs" "$targets")
		refarray+=("--nuc_min_score" "0.2")			# default Scipio score: 0.13
		refarray+=("--nuc_min_identity" "80")		# default percentage identity: 65
		refarray+=("--nuc_min_coverage" "20")		# default coverage of reference: 20
	fi
	if [ "$mode_p" = true ]; then
		refarray+=("--ptd_refs" "SeedPlantsPTD")
		refarray+=("--ptd_min_score" "0.2")			# default Scipio score: 0.2
		refarray+=("--ptd_min_identity" "75")		# default percentage identity: 65
		refarray+=("--ptd_min_coverage" "20")		# default coverage of reference: 20
	fi
	if [ "$mode_m" = true ]; then
		refarray+=("--mit_refs" "SeedPlantsMIT")
		refarray+=("--mit_min_score" "0.2")			# default Scipio score: 0.2
		refarray+=("--mit_min_identity" "75")		# default percentage identity: 65
		refarray+=("--mit_min_coverage" "20")		# default coverage of reference: 20
	fi
	if [ "$mode_e" = true ]; then
		refarray+=("--cluster_leftovers")
		refarray+=("--mmseqs_method" "easy-linclust")	# other option: "easy-cluster" for more sensitive
		refarray+=("--cl_rep_min_len" "1000")			# default min length on new ref: 500
		refarray+=("--cl_max_seq_len" "5000")			# default max length to cluster: 5000
		# note: could add more parameters if interested; see the Captus docs

		tmpdir="$(pwd)/tmp"		# temp directory for MMseqs
		# make the temp directory if it doesn't exist
		if [ ! -d "$tmpdir" ]; then
			mkdir -p "$tmpdir"
		fi
		refarray+=("--cl_tmp_dir" "$tmpdir")
	fi

	# run the extraction
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity run -H "$(pwd)" "$container" \
		extract --captus_assemblies_dir "$in_assemblies" --out "$output" "${refarray[@]}" \
 		--ram "$ram" --threads "$threads" --concurrent "$concurrent"
fi


## step 4 (alignment)
if [ "$run_step4" = true ]; then
	# defaults
	in_extractions="03_extractions"		# needs to be present
	output="./04_alignments/"
	markers="ALL"			# other options: "NUC,PTD,MIT,CLR" as comma separated list
	formats="AA,NT,GE,MA"	# AA: amino acids, NT: nucleotides, GE: genes (exon+intron),
							# GF: genes + flank, MA: matched seq, MF: matched + flank
							# note that MA and MF are applicable to CLR (clustering, mode extra)
							# also note that codon-aware alignment happens if both AA and NT
	align_method="mafft_auto"	# multiple options available; see Captus docs
	paralog_filt="naive"		# options: "naive" "informed" "both" "none"
	taper="3.0"				# default TAPER cutoff (recommended); can disable TAPER with `--disable_taper`
	clipkit="gappy"			# ClipKIT mode; other options available; see Captus docs
	clipkitg="0.9"			# ClipKIT gappyness threshold default: 0.9
	mincov="0.4"			# minimum coverage of a sequence after gaps removed default: 0.4

	# run the alignment
	srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity run -H "$(pwd)" "$container" \
		align --captus_extractions_dir "$in_extractions" --out "$output" --markers "$markers" \
		--formats "$formats" --align_method "$align_method" --filter_method "$paralog_filt" \
		--taper_cutoff "$taper" --clipkit_method "$clipkit" --clipkit_gaps "$clipkitg" \
		--min_coverage "$mincov" --ram "$ram" --threads "$threads" --concurrent "$concurrent"
fi


# record and report job completion
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)
echo -e "The Captus run finished in $duration_mins minutes or $duration_hours hours"
