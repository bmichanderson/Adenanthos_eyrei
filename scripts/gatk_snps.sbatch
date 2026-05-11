#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: May 2026
# Description: run GATK to combine GVCFs, jointly call variants, and filter SNPs to produce a VCF on Setonix
#	The input GVCFs are produced using the script gatk_hapcall_array.sbatch first
# Note: this requires a Singularity container for GATK
# Arguments:
#	-r	reference fasta file (with extension *.fasta)
#		The reference should have two associated files in the same location:
#			*.fasta.fai and *.dict
#		These will have been generated in the first step (see gatk_hapcall_array.sbatch)
#	-g	gvcf list with sampleID and /path/to/GVCF_file, tab separated and one per line
# The script will generate output in the working directory and remove the following files:
#	mydb/ *.vcf.gz *.tbi *.idx
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=gatk_snps
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --export=NONE
#SBATCH --output=gatk_snps-%j.log


# set the container
gatk_con="/software/projects/pawsey0220/anderson/singularity-containers/gatk.sif"


# load module (may need to be updated depending on Setonix updates)
module load singularity/4.1.0-slurm


# set maximum memory for Java (90% of available)
max_mem=$(echo "${SLURM_MEM_PER_NODE} * 90 / 100" | bc)m


# export the OpenMP threads variable
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"


# parse the command line
while [[ $# -gt 0 ]]		# while the number of arguments is greater than 0
do
key="$1"
case $key in
	-r)
	ref="$(readlink -f $2)"
	shift
	shift
	;;
	-g)
	gvcf_list="$(readlink -f $2)"
	shift
	shift
	;;
esac
done

# check args
if [ -z "$ref" ]; then
	echo -e "\nPlease specify a reference file with \"-r\"\n"
	exit 1
fi
if [ -z "$gvcf_list" ]; then
	echo -e "\nPlease specify a GVCF list with \"-g\"\n"
	exit 1
fi


# set up a log file
logfile="$(pwd)"/gatk_snps.log


# log start time
start="$(date +%s)"
echo -e "Starting job at $(date)\n" >> $logfile 2>&1


# generate a BED file from the fasta index, then generate "intervals"
echo -e "\nGenerating intervals\n" >> $logfile 2>&1
awk 'OFS="\t" {print $1, "0", $2}' "$ref".fai > "${ref/\.fasta/\.bed}"
singularity exec -H "$(pwd)" "$gatk_con" gatk BedToIntervalList \
	--INPUT "${ref/\.fasta/\.bed}" \
	--OUTPUT "${ref/\.fasta/\.bed}".interval_list \
	--SEQUENCE_DICTIONARY "${ref/\.fasta/\.dict}"


# run the database generation
echo -e "\nCombining GVCFs into a database\n" >> $logfile 2>&1
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$gatk_con" gatk \
	--java-options "-Xmx${max_mem} -Xms${max_mem}" GenomicsDBImport \
	--genomicsdb-workspace-path "mydb" \
	--intervals "${ref/\.fasta/\.bed}".interval_list  \
	--batch-size 16 \
	--sample-name-map "$gvcf_list" >> $logfile 2>&1


# run joint genotyping
echo -e "\nJoint genotyping...\n" >> $logfile 2>&1
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$gatk_con" gatk \
	--java-options "-Xmx${max_mem}" GenotypeGVCFs \
	--reference "$ref" \
	--variant gendb://mydb \
	--output all_variants.vcf.gz >> $logfile 2>&1


# select only SNPs for filtering
echo -e "\nSelecting SNPs from all variants...\n" >> $logfile 2>&1
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$gatk_con" gatk \
	SelectVariants \
	--variant all_variants.vcf.gz \
	--output snp_variants.vcf.gz \
	--select-type-to-include SNP >> $logfile 2>&1


# filter SNPs
echo -e "\nFiltering SNPs...\n" >> $logfile 2>&1
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$gatk_con" gatk \
	VariantFiltration \
	--reference "$ref" \
	--variant snp_variants.vcf.gz \
	--output snp_variants_filtered.vcf.gz \
	--filter-expression "QD < 6.0" --filter-name "QD6" \
	--filter-expression "SOR > 3.0" --filter-name "SOR3" \
	--filter-expression "FS > 60.0" --filter-name "FS60" \
	--filter-expression "MQ < 40.0" --filter-name "MQ40" \
	--filter-expression "MQRankSum < -12.5" --filter-name "MQRS125" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS8" >> $logfile 2>&1


# output the final VCF with only the filtered SNPs
echo -e "\nOutputting final VCF\n" >> $logfile 2>&1
srun -N 1 -n 1 -c "$SLURM_CPUS_PER_TASK" --exact singularity exec -H "$(pwd)" "$gatk_con" gatk \
	SelectVariants \
	--variant snp_variants_filtered.vcf.gz \
	--output snp_variants_final.vcf \
	--exclude-filtered >> $logfile 2>&1


# remove intermediate files
rm -r mydb/ *.vcf.gz *.tbi *.idx


# log completion
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)
echo -e "The job finished in $duration_mins minutes or $duration_hours hours\n" >> $logfile 2>&1
