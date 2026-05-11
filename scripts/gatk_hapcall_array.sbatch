#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: May 2026
# Description: run an array of GATK HaplotypeCaller jobs on Setonix
# Note: this requires two Singularity containers for Samtools and GATK
# The input is a tab-separated text file, with columns (no header):
#	sampleID	reference	bamfile
#	e.g. HERB1	/path/to/map_ref.fasta	/path/to/map_sorted.bam
# The reference fasta should have extension *.fasta!
# Note: set the "--array=__-__" variable to the number of samples and starting from 0!
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=gatk_hapcall
#SBATCH --array=0-68		# this needs to match how many samples will be analysed
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=gatk_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID} starting at $(date)"


# set the containers
samtools_con="/software/projects/pawsey0220/anderson/singularity-containers/samtools.sif"
gatk_con="/software/projects/pawsey0220/anderson/singularity-containers/gatk.sif"


# load module (may need to be updated depending on Setonix updates)
module load singularity/4.1.0-slurm


# set maximum memory for Java (90% of available)
max_mem=$(echo "${SLURM_MEM_PER_NODE} * 90 / 100" | bc)m


# export the OpenMP threads variable
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"


# check arguments
if [ -z $1 ]; then
	echo -e "Please provide a text file with sampleIDs, references and bamfiles\n"
	exit 1
fi


# Determine the sample, reference and bamfile (based on first argument text file)
samples=($(cut -f 1 $1))
refs=($(cut -f 2 $1))
bams=($(cut -f 3 $1))
sample=${samples[$SLURM_ARRAY_TASK_ID]}
ref=${refs[$SLURM_ARRAY_TASK_ID]}
ref="$(readlink -f $ref)"
bam=${bams[$SLURM_ARRAY_TASK_ID]}
bam="$(readlink -f $bam)"
echo "This job will run HaplotypeCaller for sample $sample using $SLURM_CPUS_PER_TASK CPUs"


# create a directory for the sample (if it doesn't exist) and change into it
if [ ! -d "$sample" ]; then
	mkdir -p "$sample"
fi
logfile="$(pwd)"/"$sample"/"$sample"_hapcall.log


# log start time
start="$(date +%s)"
echo -e "Starting job at $(date)\n" >> $logfile 2>&1


# generate an index for the reference if it doesn't exist
if [ ! -f "${ref/\.fasta/\.dict}" ]; then
	singularity exec -H "$(pwd)" "$gatk_con" gatk CreateSequenceDictionary --REFERENCE $ref
fi
if [ ! -f "$ref".fai ]; then
	singularity exec -H "$(pwd)" "$samtools_con" samtools faidx $ref
fi


# change to the sample directory
cd "$sample"


# specify read groups
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$gatk_con" \
	gatk AddOrReplaceReadGroups \
	--INPUT "$bam" \
	--OUTPUT "$sample"_rg.bam \
	--RGLB "mylib" \
	--RGPL ILLUMINA \
	--RGPU "myunit" \
	--RGSM "$sample" \
	--CREATE_INDEX true \
	>> $logfile 2>&1


# run the calling
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$gatk_con" \
	gatk --java-options "-Xmx${max_mem}" HaplotypeCaller \
	--reference "$ref" \
	--input "$sample"_rg.bam \
	--output "$sample".gvcf.gz \
	--emit-ref-confidence GVCF \
	>> $logfile 2>&1


# log completion
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
echo -e "The job finished in $duration_mins minutes\n" >> $logfile 2>&1
cd ..


# Echo the end of the current array job
echo "Job ${SLURM_JOB_ID} for sample $sample finished at $(date)"
