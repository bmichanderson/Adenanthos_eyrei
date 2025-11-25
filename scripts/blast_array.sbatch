#!/bin/bash --login

######################
# Author: B.M. Anderson
# Date: Aug 2024
# Modified: Apr 2025 (output in the current directory rather than sample directories; no longer use sample; change BLAST output)
# Description: run an array of BLAST jobs on Setonix, using a text file (arg1) with queries and subjects
# Note: the text file should provide paths to query and subject fasta file(s) (tab separated, one pair per line)
#	i.e.	/path/to/query1 \t /path/to/subject1
#			/path/to/query2 \t /path/to/subject2
# Note: the BLAST search will only keep hits longer than 500 bp and >90 identity (can be changed below)
# Note: it assumes the queries are unique and will name output using everything before the first underscore in the query names
# Note: set the "--array=__-__" variable to the number of queries and starting from 0!
######################


# Assign slurm variables
#SBATCH --account=pawsey0220
#SBATCH --job-name=blast
#SBATCH --array=0-9		# this needs to match the number of blast calls
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --output=blast_array-%j.log


# Echo the current array job
echo "This is job ${SLURM_JOB_ID} with task ID ${SLURM_ARRAY_TASK_ID}"


# Load the appropriate module (may change with updates)
module load blast/2.12.0--pl5262h3289130_0


# check arguments
if [ -z $1 ]; then
	echo -e "Please provide a text file with paths to query and subject fastas\n"
	exit 1
fi


# Determine the query and subject files (based on first argument text file)
queries=($(cut -f 1 $1))
query=${queries[$SLURM_ARRAY_TASK_ID]}
subjects=($(cut -f 2 $1))
subject=${subjects[$SLURM_ARRAY_TASK_ID]}


# Grab a prefix from the name of the query (assumes there is an underscore)
filename=$(basename $query)
prefix=${filename/_*/}


# Run BLAST for the given input
echo "This job is for query $(basename $query) and subject $(basename $subject)"

makeblastdb -dbtype nucl -in "$subject" -out dbnucl_"$SLURM_ARRAY_TASK_ID" > /dev/null

echo "qacc qlen sacc slen pident length qstart qend sstart send" | tr ' ' '\t' > "$prefix"_blast_out.tab

blastn -query "$query" -db dbnucl_"$SLURM_ARRAY_TASK_ID" \
	-outfmt "6 qacc qlen sacc slen pident length qstart qend sstart send" \
	-max_target_seqs 10 | awk '$5 > 90' | awk '$6 > 500' >> "$prefix"_blast_out.tab

rm dbnucl_"$SLURM_ARRAY_TASK_ID".*

echo "The job for query $(basename $query) finished"
