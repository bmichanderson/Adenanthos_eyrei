#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Jun 2024
# Modified: Mar 2025 (made deduplication optional)
# Description: run a single Illumina QC task on Setonix
# Note: this should be called by the array script
######################


# set max memory to use (85% of total for Java)
max_mem=$(echo "${SLURM_MEM_PER_NODE}*85/100" | bc)m
echo "This job is going to use $SLURM_CPUS_PER_TASK cores with max mem of $max_mem"


# define software location
bbmap_loc="/software/projects/pawsey0220/anderson/software/bbmap"


# set whether to run deduplication or not
dedup="True"


# bbduk args
min_len="50"		# minimum length of reads to keep after trimming
adapter_r1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_r2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
kmer="21"		# kmer value to use when searching for adapter in entire read
mismatch="2"		# how different an adapter sequence kmer can be for the kmer match
minq="20"		# minimum quality to trim to
kshort="8"
adapter_short1="${adapter_r1:0:15}"
adapter_short2="${adapter_r2:0:15}"
forceleft="9"	# trim all reads this many bases at the start of the read; if weirdness detected


# parse the command line
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
esac
done


# check args
if [ -z "$read1" ] || [ -z "$read2" ]; then		# if input files are not specified
	exit 1
fi


## start QC
# set up a folder and log file named based on the read file name
file_name=$(basename ${read1})
prefix=${file_name/_*/}
if [ ! -d "$prefix" ]; then
	mkdir -p "$prefix"
fi
cd "$prefix"
logfile="$prefix".log


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting Illumina QC at $(date)\n**********" >> $logfile 2>&1


# run readlength stats on the raw files
echo -e "\n\n**********\nAssessing raw readlength distributions\n**********\n" >> $logfile 2>&1
"$bbmap_loc"/readlength.sh in1="$read1" in2="$read2" out=readlength_raw.txt -Xmx"$max_mem" >> $logfile 2>&1


# remove duplicates (optical; dupedist is for NovSeq), if requested
if [ "$dedup" == "True" ]; then
	echo -e "\n\n**********\nRemoving optical duplicates with clumpify\n**********\n"  >> $logfile 2>&1
	"$bbmap_loc"/clumpify.sh in1="$read1" in2="$read2" out1=temp_1.fastq.gz out2=temp_2.fastq.gz \
		dedupe dupedist=12000 optical subs=0 -Xmx"$max_mem" >> $logfile 2>&1
else
	echo -e "\n\n**********\nSkipping optical deduplication\n**********\n"  >> $logfile 2>&1
	ln -s "$read1" ./temp_1.fastq.gz
	ln -s "$read2" ./temp_2.fastq.gz
fi


# trim adapters, for quality, and PhiX
echo -e "\n\n**********\nTrimming adapters, for quality, and PhiX\n**********\n" >> $logfile 2>&1
"$bbmap_loc"/bbduk.sh in1=temp_1.fastq.gz in2=temp_2.fastq.gz out1=temp1_1.fastq.gz out2=temp1_2.fastq.gz ktrim=r \
	literal="$adapter_r1","$adapter_r2" k="$kmer" hdist="$mismatch" minlength="$min_len" mink=15 \
	tbo tpe forcetrimleft="$forceleft" -Xmx"$max_mem" >> $logfile 2>&1
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/bbduk.sh in1=temp1_1.fastq.gz in2=temp1_2.fastq.gz out1=temp2_1.fastq.gz out2=temp2_2.fastq.gz ktrim=r \
	literal="$adapter_short1","$adapter_short2" k="$kshort" restrictright="$kmer" minlength="$min_len" \
	qtrim=r trimq="$minq" -Xmx"$max_mem" >> $logfile 2>&1
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/bbduk.sh in1=temp2_1.fastq.gz in2=temp2_2.fastq.gz out1=trim_1.fastq.gz out2=trim_2.fastq.gz \
	ref=phix k=31 hdist="$mismatch" -Xmx"$max_mem" >> $logfile 2>&1
rm temp*


# correct sequencing errors
echo -e "\n\n**********\nCorrecting sequencing errors\n**********\n" >> $logfile 2>&1
"$bbmap_loc"/reformat.sh in1=trim_1.fastq.gz in2=trim_2.fastq.gz out=temp.fastq.gz \
	-Xmx"$max_mem" >> $logfile 2>&1
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/bbmerge.sh in=temp.fastq.gz out=ecco.fastq.gz ecco mix vstrict ordered \
	prefilter=2 prealloc=t minlength="$min_len" -Xmx"$max_mem" >> $logfile 2>&1
mv ecco.fastq.gz temp.fastq.gz
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/clumpify.sh in=temp.fastq.gz out=eccc.fastq.gz ecc passes=4 -Xmx"$max_mem" >> $logfile 2>&1
mv eccc.fastq.gz temp.fastq.gz
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/tadpole.sh in=temp.fastq.gz out=ecct.fastq.gz ecc k=62 ordered \
	prefilter=2 prealloc=t errormult1=64 threads="$SLURM_CPUS_PER_TASK" -Xmx"$max_mem" >> $logfile 2>&1
mv ecct.fastq.gz temp.fastq.gz
echo -e "\n\n" >> $logfile 2>&1
"$bbmap_loc"/bbduk.sh in=temp.fastq.gz out1=correct_1.fastq.gz out2=correct_2.fastq.gz \
	minlength="$min_len" -Xmx"$max_mem" >> $logfile 2>&1
rm temp* trim*


# rename the output
mv correct_1.fastq.gz "$prefix"_R1.fastq.gz
mv correct_2.fastq.gz "$prefix"_R2.fastq.gz


# run readlength stats on the final files
echo -e "\n\n**********\nAssessing final readlength distributions\n**********\n" >> $logfile 2>&1
"$bbmap_loc"/readlength.sh in1="$prefix"_R1.fastq.gz in2="$prefix"_R2.fastq.gz out=readlength_qc.txt \
	-Xmx"$max_mem" >> $logfile 2>&1
echo -e "\n\n" >> $logfile 2>&1


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)
echo -e "\nFinished Illumina QC at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours\n" >> $logfile 2>&1
echo "Job ${SLURM_ARRAY_TASK_ID} for sample "$prefix" finished in $duration_mins minutes"
