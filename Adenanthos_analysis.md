Genome skimming for *Adenanthos*
================================
B.M. Anderson
started October 2023
new data available August 2024  
new data available March 2025  

These notes describe the assembly and analysis of *Adenanthos* genome skimming data generated to sort out the taxonomic status of *Adenanthos eyrei* E.C.Nelson  

# Data
The only known collection of *Adenanthos eyrei* (the holotype E.C. Nelson ANU 17044 -- CANB 237043) was sampled  
DNA extractions used a modified CTAB method and were pooled  
Library preparation used the xGen cfDNA & FFPE DNA Library Prep v2 MC Kit (10010207) from Integrated DNA Technologies  
Sequencing was done on an Illumina NovaSeq X Plus in the 150 bp paired-end format  
Raw reads were obtained as two pairs of two read files (four files)  
A putative hybrid sample (E.C. Nelson ANU 16981 -- PERTH 01711202) was later also sequenced the same way  

Further sequencing of fresh samples was done similarly, but with less effort (smaller files)  
Library preparation used the Illumina DNA Prep (M) Tagmentation kit (20060059) (using Nextera adapters)  

All data is now stored on the DBCA Microsoft Azure Storage account "floragenetics"  
To view file info, use the Azure Command Line Interface
```bash
az login --use-device-code
# use a browser to go to https://microsoft.com/devicelogin and enter the code you are given in the terminal; hit enter after for default
# check contents
az storage blob list \
--account-name floragenetics \
--container-name adenanthos-eyrei \
--output table \
--auth-mode login > data.txt
```

Cut the filenames from the resulting text file
```bash
grep fastq.gz data.txt | cut -f 1 -d " " | cut -f 2 -d "/" > files.txt
```
The list of file names can be uploaded to Setonix for transfer  

Following advice, I have downloaded the AzCopy binary to my software folder on Setonix: `/software/projects/pawsey0220/anderson/software/azcopy/azcopy`, then I created a symbolic link to it in my `~/bin/` directory, so I can launch anywhere with `azcopy`  

Create a `samples.txt` file with the names of all samples (should match read file names before first underscore), one per line  
Upload it to Setonix under a directory `Adenanthos`  

To copy the files in `files.txt`, run an interactive job from a new folder `Adenanthos/raw` with the `files.txt` file in it
```bash
salloc -p copy -n 1 -N 1 -c 4 --mem=8G -A pawsey0220 --time=01:00:00
# wait for it to start...
keyctl session workaround		# starts a new keyctl session
azcopy login
# use a browser to go to https://microsoft.com/devicelogin and enter the code you are given in the terminal
azcopy copy "https://floragenetics.blob.core.windows.net/adenanthos-eyrei/rawdata-AGRF-shotgun/*" . --list-of-files files.txt

# after copy is finished, concatenate appropriate files together, using the samples.txt file
# first, change the name of the files with prefix "Adenanthos" to be consistent
for file in Adenanthos*.fastq.gz; do mv $file ${file/Adenanthos/}; done
for sample in $(cat ../samples.txt); do
cat "$sample"_*R1.fastq.gz > "$sample"_R1_comb.fastq.gz && rm "$sample"_*R1.fastq.gz
cat "$sample"_*R2.fastq.gz > "$sample"_R2_comb.fastq.gz && rm "$sample"_*R2.fastq.gz
mv "$sample"_R1_comb.fastq.gz "$sample"_R1.fastq.gz
mv "$sample"_R2_comb.fastq.gz "$sample"_R2.fastq.gz
done

# end the session and the interactive job
exit
exit
```

## Quality assessment and trimming
Run FastQC v. 0.11.9 on the raw data  
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/fastqc.sbatch *.gz
```

Results are different for two groups of samples:  
The herbarium samples (HERB1 and HERB2) show high amounts of Illumina adapter, non-random bases for the first 9 bp, and some issues with quality at the ends of reads  
The fresh samples (the rest) typically show small amounts of Nextera Transposase adapter and non-random bases for the first 14 bp  

Run Illumina QC steps, removing the first 9/14 bp of reads, trimming adapters and correcting errors  
(from a new `Adenanthos/qc` folder)
```bash
ln -s /scratch/pawsey0220/anderson/Adenanthos/raw/*.gz .

# for the fresh samples (74 samples), run an array of jobs
rm HERB*.gz
cp /software/projects/pawsey0220/anderson/scripts/illumina_qc.sh ./illumina_qc_fresh.sh
# modify the copy of illumina_qc.sh to trim 14 bp, turn off deduplication, and include Nextera adapters:
# adapter_r1="CTGTCTCTTATACACATCTGACGCTGCCGACGA"
# adapter_r2="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
cp /software/projects/pawsey0220/anderson/scripts/illumina_qc_array.sbatch .
# modify the array sbatch to match how many samples (pairs of files = 74, so set 0-73) to analyse
# also update it to execute the script illumina_qc_fresh.sh
sbatch illumina_qc_array.sbatch

# for the two herbarium samples, run separately afterwards
ln -s /scratch/pawsey0220/anderson/Adenanthos/raw/HERB*.gz .
cp /software/projects/pawsey0220/anderson/scripts/illumina_qc.sh ./illumina_qc_HERB.sh
# modify the script to turn off deduplication
cp illumina_qc_array.sbatch illumina_qc_array_HERB.sbatch
# modify the array sbatch to match how many samples (pairs of files = 2, so set 0-1) to analyse
# also update it to execute the script illumina_qc_HERB.sh, use 40 cores and only look for HERB*.gz files
sbatch illumina_qc_array_HERB.sbatch
```

Summarise the results of cleaning and filtering  
```bash
echo -e "Sample\tReads_raw\tRemoved\tReads_qc" > qc_summary.txt
for prefix in $(cat ../samples.txt); do
paste <(echo "$prefix") <(grep "Reads" "$prefix"/readlength_raw.txt | cut -f 2) \
<(grep "Total Removed" "$prefix"/"$prefix".log | cut -f 2 | cut -f 1 -d " " | paste -sd+ | bc) \
<(grep "Reads" "$prefix"/readlength_qc.txt | cut -f 2) >> qc_summary.txt
done
```

Run another FastQC
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/fastqc.sbatch */*.gz
```

## Downsampling
Given the large number of reads present for the herbarium specimens (c. 170 M read pairs), downsample for the initial SPAdes assembly  
To do so, specify an approximate number of reads (80 M = 40 M read pairs) and prioritize longer reads  
```bash
for sample in HERB1 HERB2; do
cd "$sample"/
sbatch /software/projects/pawsey0220/anderson/scripts/downsample.sbatch -1 "$sample"_R1.fastq.gz \
-2 "$sample"_R2.fastq.gz -a 80000000 -r yes
cd ..
done
```
This resulted in c. 73 M reads (average length 140 bp) for HERB1 and c. 79 M reads (average length 119 bp) for HERB2  
HERB2 was apparently more fragmented, and it has a greater number of short reads  

# Overall Assembly
Use SPAdes v. 4.0.0 to run assemblies for all samples  
These can be used as input to GetOrganelle or to extract contigs that match high copy targets like the nuclear ribosomal region  

From a directory `Adenenthos/assembly` with links to the cleaned reads, run an assembly script in an array  
```bash
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz .
rm HERB{1,2}_R{1,2}.fastq.gz	# remove links to the full size files, leaving the downsampled ones
cp /software/projects/pawsey0220/anderson/scripts/spades.sh .
cp /software/projects/pawsey0220/anderson/scripts/spades_array.sbatch .
# modify the array script to match how many samples to assemble (76, so 0-75)
sbatch spades_array.sbatch
```
Took from 69 to 205 min per sample using 32 cores  

If there are any errors, can re-try (more memory) with an interactive job
```bash
salloc -p work -n 1 -N 1 -c 32 --mem=128G --time=06:00:00 -A pawsey0220
./spades.sh -1 GLA1-04*R1.fastq.gz -2 GLA1-04*R2.fastq.gz -o GLA1-04
exit
```
This sample took about 4 hours  

# Plastome
Assemble full plastomes for all samples, then use alignments to assess similarity  

## Assembly
Assemble a first draft of each plastome using the following approach:  
1) use GetOrganelle v. 1.7.7.1 to pull out and clean putative plastome assemblies from the SPAdes assemblies  
2) combine all recovered circular assemblies into a single reference file  
3) map reads to the reference for each sample, then assemble mapped reads with Unicycler v. 0.5.1  
4) examine Unicycler assemblies in Bandage v. 0.8.1 to detect and extract expected tripartite structure  
5) re-orient all plastomes to start from the LSC  
6) map reads for each sample to the respective assembly to detect problems  
7) correct problems and annotate final assemblies using Chloe  

### GetOrganelle
From a new directory `Adenanthos/plastome`, download the GetOrganelle databases in preparation for assembly  
```bash
salloc -p copy -n 1 -N 1 -c 4 --mem=8G -A pawsey0220 --time=01:00:00
module load singularity/4.1.0-slurm
container="/software/projects/pawsey0220/anderson/singularity-containers/getorganelle.sif"
singularity exec -H "$(pwd)" "$container" get_organelle_config.py --config-dir "$(pwd)"/.GetOrganelle --add embplant_pt,embplant_mt
# wait until it finishes (about 30 sec)
exit
```
This puts the databases in the folder `.GetOrganelle` in the current directory  

Create a list of folder paths to the SPAdes assemblies, then launch plastome assemblies  
```bash
echo /scratch/pawsey0220/anderson/Adenanthos/assembly/*/ | tr -s ' ' '\n' > spades_assemblies.txt
cp /software/projects/pawsey0220/anderson/scripts/plastome_from_assembly.sh .
cp /software/projects/pawsey0220/anderson/scripts/plastome_array.sbatch .
# modify the array script for the number of samples (76, so specify 0-75)
sbatch plastome_array.sbatch spades_assemblies.txt
```
Took about 4 to 40 (84) minutes per sample using 16 cores  

Check how many recovered a circular genome, how many contigs/paths, and lengths
```bash
for sample in $(cat ../samples.txt); do
paste <(echo $sample) \
<(grep "Result status" ${sample}/getorganelle_out/get_org.log.txt | cut -f 5 -d ":" | sed 's/^[[:space:]]*//') \
<(ls ${sample}/getorganelle_out/*.fasta | wc -l) \
<(grep "Total:LSC:SSC:Repeat" ${sample}/getorganelle_out/get_org.log.txt | cut -f 2 -d "=" | sed 's/^[[:space:]]*//') >> summary.txt
done
```
Most (66/76) of the samples returned a complete circle  

### Reference generation
Use the circular fastas from GetOrganelle to create a reference for mapping and assembly with Unicycler  
The fasta reference names are long, so change them to simple numbers  
Compile them into a single initial file
```bash
ls */getorganelle_out/*complete.graph1.1*.fasta > fasta_list.txt
index=1
for fasta in $(cat fasta_list.txt); do
cat $fasta | sed "s/^>.*$/>$index/g" >> reference.fasta
index=$((index + 1))
done
```

Download the `reference.fasta` to reformat it using matches to known plastid genes  
Pull out the first fasta and annotate with Chloe (either via Singularity container or at GeSeq online)  
```bash
grep ">" -A1 reference.fasta | head -n 2 > 1.fasta
singularity run -H "$(pwd)" ~/singularity-containers/chloe.sif annotate --reference=/chloe_references -r cp --gbk --no-gff 1.fasta
```
Extract two reference gene sequences from `1.chloe.gbk` for orienting consistently  
(in order for the script to work, `  ORGANISM  Adenanthos` needs to be manually added as a second line to the Chloe output, and the first line needs to have spaces adjusted -- see Annotation section below)
*psbA* is typically near the start of the LSC and in reverse orientation  
*ccsA* is typically in the SSC and in the forward orientation  
```bash
python3 ~/scripts/genbank_parse.py -g psba 1.chloe.gbk
mv Adenanthos_extract.fasta psba.fasta
python3 ~/scripts/genbank_parse.py -g ccsa 1.chloe.gbk
mv Adenanthos_extract.fasta ccsa.fasta
```

To enable alignment and consensus generation, orient assemblies consistently  
The script to do this requires BLAST installed  
```bash
python3 ~/scripts/fasta_splitting.py reference.fasta
for ref in {1..9}*.fasta; do
bash ~/scripts/reform_cp.sh "$ref" psba.fasta ccsa.fasta
mv new_cp.fasta "$ref"
done
```

Concatenate together, then run an alignment with MAFFT  
After that finishes, generate a consensus sequence of the alignment
```bash
cat {1..9}*.fasta > reference.fasta
mafft --auto --thread 8 reference.fasta > aligned.fasta
# manually check that it looks OK
python3 ~/scripts/consensus.py -t 0.1 aligned.fasta
mv consensus.fasta reference.fasta && rm aligned.fasta {1..9}*.fasta
```
Upload the new `reference.fasta` to the folder with plastome assemblies (`Adenanthos/plastome`)

### Unicycler
Ideally, Unicycler is an effective tool for assembling circular structures using the combination of long and short reads  
However, even with only short reads (as for this dataset) Unicycler can optimise SPAdes assemblies to pull out circular structures  

For each sample, map its reads against the reference fasta, then run Unicylcer to assemble the mapped reads and refine the circles  
For the mapping, the Unicycler script will adhere to semi-strict settings for BBMap (see the script for details)  
It has the option to normalise coverage to 1000x after mapping prior to SPAdes assemblies, but keep that off for now  
Link the QC reads into the `Adenanthos/plastome` folder, then copy over and run Unicycler scripts
```bash
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz .
# remove the full sized HERB files (use the downsampled set of available reads)  
rm HERB{1,2}_R{1,2}.fastq.gz
# create the file with sample names and reference path
for sample in $(cat ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/reference.fasta) >> unicycler_file.txt
done
# copy over the scripts
cp /software/projects/pawsey0220/anderson/scripts/unicycler.sh .
cp /software/projects/pawsey0220/anderson/scripts/unicycler_array.sbatch .
# modify the array script to match how many samples to assemble (76, so 0-75)
sbatch unicycler_array.sbatch unicycler_file.txt
```
The jobs took about 5–21 minutes each  

Download the `assembly.gfa` files for each sample  
```bash
for sample in $(cat ../samples.txt); do
cp $sample/unicycler/assembly.gfa ./"$sample"_assembly.gfa
done
# then download the gfa files
```

Use Bandage to examine the assemblies (*.gfa files) and extract circular structures  
The typical tripartite structure will appear as two circular regions (LSC and SSC) connected by a single contig (the IR region)  
Using display of both strands, Bandage can be used to extract a path that traverses the regions  
It should start on one of the LSC strands and cross the IR twice  
If there are multiple paths available, choose the higher depth option  
Output as `{sample}_assembly.fasta`  

Some samples showed indications of lower depth alternative sequences (deleted to form the circles), which may indicate contamination or transfer to another part of the cell  

A complete circle could not be recovered for one sample (FOR3-01)  
Run again, this time mapping to the circular assembly recovered by GetOrganelle for that sample (rather than the consensus reference)  
The second run recovered a complete circle  

### Mapping
Map reads for all samples to assess how clean the assemblies are and if there are indications of problems  

First, reform each of the assemblies using the references genes
```bash
for sample in $(cat ../samples.txt); do
bash ~/scripts/reform_cp.sh "$sample"_assembly.fasta psba.fasta ccsa.fasta
sed "s/^>.*/>${sample} plastome/" new_cp.fasta > "$sample"_assembly.fasta
rm new_cp.fasta
done
```

For mapping and ensuring there is coverage across the boundaries between the IR and the LSC, add 300 bp to each end  
Create a copy of the plastome (remembering that any coordinates will need to be adjusted accordingly)  
```bash
for sample in $(cat ../samples.txt); do
python3 ~/scripts/extend_circle.py "$sample"_assembly.fasta -b 300
mv new_contig.fasta "$sample"_assembly_map.fasta
done
```

Upload to a new folder on Setonix `Adenanthos/plastome/map`  
In the new folder, create links to the QC reads again, removing the non-downsampled versions of the HERB samples
```bash
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz .
rm HERB{1,2}_R{1,2}.fastq.gz
```

Make a text file with each sample name followed by reference fasta file path
```bash
for sample in $(cat ../../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/${sample}_assembly_map.fasta) >> mapping_file.txt
done
```

Copy over the scripts, adjust for the correct number of samples in the array script, then launch mapping
```bash
cp /software/projects/pawsey0220/anderson/scripts/mapping.sh .
cp /software/projects/pawsey0220/anderson/scripts/mapping_array.sbatch .
# modify the array script to match how many samples to map (76, so 0-75)
sbatch mapping_array.sbatch mapping_file.txt
```
Took about 20 sec to 3 min per sample  

Run a summary of the mapping using Samtools
```bash
cp /software/projects/pawsey0220/anderson/scripts/samtools_summary.sh .
cp /software/projects/pawsey0220/anderson/scripts/samtools_array.sbatch .
# modify the array script to match how many samples to map (76, so 0-75)
# ensure the samtools_summary.sh uses min_depth=15, het_fraction=0.4, call_fraction=0.6
sbatch samtools_array.sbatch -f ../../samples.txt -e 300
```
Took about 3 sec to 1 min per sample  

Concatenate the summary files together for input to a spreadsheet
```bash
index=1
for sample in $(cat ../../samples.txt); do
if [ $index == 1 ]; then
	head -n 1 "$sample"/"$sample"_samtools_summary.txt > samtools_overall.txt
fi
tail -n 1 "$sample"/"$sample"_samtools_summary.txt >> samtools_overall.txt
index=$((index + 1))
done
```

Download the mapping files (*.bam and *.bai and the reference fastas) to visually inspect sites of suspected hets or problems  
Use the Integrated Genomics Viewer or similar program for viewing BAM files and mapping references  

Most assemblies appear clean, but a few show indications of a secondary (lower depth) sequence near 103–109 and 136–142 kbp  
This corresponds to a region of the IR  
Run an alignment of consensus sequences
```bash
cat */*consensus.fasta > input.fasta
mafft --auto --thread 8 input.fasta > aligned.fasta
``` 
Examining the alignment indicates the het sites are not phylogenetically informative for the study group  
(i.e. no match or variation at that site besides the hets in samples with them)  
Looking at other samples without those hets called suggests the sequence is sometimes present at lower depth  
It is possible the sequence represents a transfer to the mitogenome or nuclear genome  

One site in the SSC consistently appeared to be heterozygous in many samples, so it was scored as heterozygous even when found below the threshold (0.4) in other samples  

One issue concerns populations SER1 and SER2, where there are indications of an assembly error near the junction of the SSC and IRs  
(using the extended mapping coordinates, so subtract 300 for the assembly coordinates)  

113160–113190  
SER1-01		TGTTTCATTTCTTCTGGAAAAATCACAAAGA  
(read)		TGTTTCATTTCTTCTGGAAAAATCA-AAAGA  
(read)		                        A-AAAGAAATTCAAAATTAAATATAAATAAAAAAGTAGTAAATTACTAAAAAGATTAATACATAAGATAAATA  
(read)		TGTTTCATTTCTTCTGGAAAAATCACAAACACTTTTTTGATT  

132540–132560  
SER1-01				 TCAAAAAAGTGTTT-TGATTTT  
(read)				 TCAAAAAAGTGTTTGTGATTTT  
(read)	TATTTATATTTAATTTTGAATTTCTTT-TGATTTT  

So, need to remove a C and insert a G as follows:  
SER1-01, remove C at 113185 and insert G after 132553 (T)  
SER1-02, remove C at 113185 and insert G after 132553 (T)  
SER1-03, remove C at 113185 and insert G after 132553 (T)  
SER1-04, remove C at 113210 and insert G after 132576 (T)  

SER2-01, remove C at 113157 and insert G after 132523 (T)  
SER2-02, remove C at 113157 and insert G after 132523 (T)  
SER2-03, remove C at 113178 and insert G after 132544 (T)  
SER2-04, remove C at 113178 and insert G after 132544 (T)  

Additionally, SER1-05 and SER2-05 had heterozygous sites dispersed throughout the plastome,  
suggesting possible contamination at lower depth  

To make adjustments (het site in all samples, assembly corrections), create a text file from a spreadsheet with columns (no header):  
fasta_id	position(1-based)	current_base	new_base  

(remember to adjust by 300!)  

Since the assembly fasta files are named as: ">sample plastome", use sample as the first column in the text file (`positions.txt`)  
To delete a base, use 'C' '-'; to insert, use 'T' 'TG'  
```bash
for fasta in *_assembly.fasta; do
python3 ~/scripts/bp_alter.py -f $fasta -p positions.txt
mv output_alter.fasta ${fasta/assembly/mod}
done
```

If wanting to check the new assemblies, create new mapping references
```bash
for sample in $(cat ../samples.txt); do
python3 ~/scripts/extend_circle.py "$sample"_mod.fasta -b 300
mv new_contig.fasta "$sample"_mod_map.fasta
done
```

Upload and repeat the mapping as before  

## Annotation
Options to annotate include using the progrma Chloe (https://github.com/ian-small/chloe) in a Singularity container,  
or to run via the web interface at GeSeq (https://chlorobox.mpimp-golm.mpg.de/geseq.html)  

If using the website, prep the fasta files for upload
```bash
for file in *assembly.fasta; do sed 's/ plastome//' $file > ${file/_assembly/}; done
```
When uploading to GeSeq, select to support all annotations with Chloe v0.1.0 ("CDS" "tRNA" "rRNA"); turn off BLAT settings  
The downloaded job has individual GenBank files included (copy them from the downloaded folder, and remove the job name from file names)  
```bash
for sample in $(cat ../samples.txt); do mv *"$sample"*.gb "$sample".gb; done
```

*Preferred option*
If using Chloe via the Singularity container (latest development version -- unversioned), the output needs to be adjusted after running  
Chloe can output GenBank and EMBL flat files using the GenomicAnnotations Julia package (https://github.com/BioJulia/GenomicAnnotations.jl)  
First, run the annotation, generating GenBank and EMBL files
```bash
singularity run -H "$(pwd)" ~/singularity-containers/chloe.sif annotate \
--reference=/chloe_references --no-transform -r cp --gbk --embl --no-gff *assembly.fasta
```

Adjust the header lines in the GenBank files to the format expected by Biopython (**NOTE: this is no longer required if Ian accepts my pull request**)
```bash
IFS=" "
for gbfile in *.gbk; do
myarray=( $(head -n 1 $gbfile | tr -s ' ') )
len_name=$(printf ${myarray[1]} | wc --chars)
len_seq=$(printf ${myarray[2]} | wc --chars)
spaces=$(echo "28 - $len_name - $len_seq" | bc)
out_spaces=$(printf '%*s' $spaces)
len_type=$(printf ${myarray[5]} | wc --chars)
spaces2=$(echo "9 - $len_type" | bc)
out_spaces2=$(printf '%*s' $spaces2)
new_line=$(echo "${myarray[0]}       ${myarray[1]}${out_spaces}${myarray[2]} ${myarray[3]}     ${myarray[4]}    ${myarray[5]}${out_spaces2}${myarray[6]} ${myarray[7]}")
sed -i "s/LOCUS.*$/$new_line/" $gbfile
done
IFS=$' \t\n'
```
The GenBank flat files will now be parseable  

Add organism info to the flat files for submission and parsing if needed
```bash
ls *assembly.fasta | cut -f1 -d "_" > samples.txt
for sample in $(grep CUN samples.txt); do
sed -i '1 a SOURCE      chloroplast Adenanthos cuneatus\n  ORGANISM  Adenanthos cuneatus' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos cuneatus"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep DOB samples.txt); do
sed -i '1 a SOURCE      chloroplast Adenanthos dobsonii\n  ORGANISM  Adenanthos dobsonii' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos dobsonii"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep FOR samples.txt); do
sed -i '1 a SOURCE      chloroplast Adenanthos forrestii\n  ORGANISM  Adenanthos forrestii' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos forrestii"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep GLA samples.txt); do
sed -i '1 a SOURCE      chloroplast Adenanthos glabrescens\n  ORGANISM  Adenanthos glabrescens subsp. exasperatus' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos glabrescens subsp. exasperatus"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep HERB samples.txt); do
sed -i '1 a SOURCE      chloroplast Adenanthos eyrei\n  ORGANISM  Adenanthos eyrei' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos eyrei"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep ILE samples.txt); do
sed -i '1 a SOURCE      chloroplast Adenanthos ileticos\n  ORGANISM  Adenanthos ileticos' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos ileticos"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep SER samples.txt); do
sed -i '1 a SOURCE      chloroplast Adenanthos sericeus\n  ORGANISM  Adenanthos sericeus' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos sericeus"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
```

Replace the unknown lengths with the number at the end of the embl files (could be avoided with an improvement to Chloe code perhaps)
```bash
for file in *.embl; do
mynum=$(grep Sequence $file | tr -s ' ' | cut -f 3 -d ' ')
sed -i "s/??/${mynum}/" $file
done
```

Submitting any plastomes to GenBank requires altering the GenBank flat files by extracting the features into a table with specific formatting  
This is potentially a nuisance (tools are available to convert these files, e.g. https://chlorobox.mpimp-golm.mpg.de/GenBank2Sequin.html)  
There is the additional difficulty that Chloe omits the required feature qualifier "product", and implements "number" of introns differently  
**Note: this is now at least partly corrected with my next pull request (if Ian accepts it)**  
It also adds extraneous ones ("name", "ID", "parent") that are applicable to the GFF3 format, and adds "locus_tag", which for GenBank is added later  
An alternative route is to submit to ENA instead of GenBank, since ENA accepts EMBL flat files (similar format to GenBank)  
ENA asks that a locus_tag be provided, and one needs to be registered with the project before submitting (this could be search replaced perhaps)  

Remove extraneous lines from the GenBank and ENA files
```bash
for sample in $(cat samples.txt); do
sed -i '/ID=/d' "$sample"_*.gbk
sed -i '/ID=/d' "$sample"_*.embl
sed -i '/Parent/d' "$sample"_*.gbk
sed -i '/Parent/d' "$sample"_*.embl
sed -i '/Name/d' "$sample"_*.gbk
sed -i '/Name/d' "$sample"_*.embl
done
```

## Alignments
One option is to align across the entire plastome, another is to use only the CDS regions from the GenBank flat files  

### Full plastomes
Align across all plastomes (use the fasta assembly files that were used for annotation)
```bash
cat *assembly.fasta > input.fasta
mafft --thread 8 --auto input.fasta > aligned.fasta
```

### CDS
Extract CDS from the GenBank flat files
```bash
# list CDS in all GenBank files
for file in *.gbk; do python3 ~/scripts/genbank_parse.py -l CDS $file >> temp; done
sort temp | uniq > cds_list.txt
rm temp
# extract the nucleotide alignments of the CDS
for file in *.gbk; do
python3 ~/scripts/genbank_parse.py -f cds_list.txt $file &&
mv Adenanthos*_nucl_extract.fasta ${file/\.chloe\.gbk/_nucl_extract\.fasta}
done
# check duplicated (or more) regions (expected for ndhB, rpl2, rpl23, rps7, rps12, ycf2 in the IR)
grep ">*_[2,3,4]" *assembly_nucl_extract.fasta | cut -f 2 -d ">" | cut -f 1 -d " " | sort | uniq > dup_list.txt
```

Remove all second copies (with "_2")
```bash
python3 ~/scripts/remove_fastas.py _2 *nucl_extract.fasta
for file in mod*.fasta; do mv "$file" "${file/mod_/}"; done
```

Compile the extracts into individual gene files, prior to aligning them
```bash
# run this in a separate folder "CDS_alignment"
python3 ~/scripts/extract_sort_comb.py ../*extract.fasta
```
This produced 79 gene files, each with 76 samples (no missing)  

Prior to continuing with alignment, rename the entries to just the sample name  
```bash
# from the "CDS_alignment" folder
for file in *.fasta; do
gene="${file/.fasta/}"
sed "/>/s/${gene} from //g" $file | sed "/>/s/ .*$//" > temp
mv temp $file
done
```

Run translation, alignment, then back substitute the nucleotides with `pal2nal.pl` (in Singularity container)
```bash
python3 ~/scripts/translate.py -c 11 *.fasta
# remove stop codons at the ends of genes
sed -i 's/\*$//g' *_prot.fasta
# align
for translation in *_prot.fasta; do
mafft --auto --thread 8 "$translation" > "${translation/.fasta/.aln}"
done
# convert back
for alignment in *_prot.aln; do
singularity exec -H "$(pwd)" ~/singularity-containers/phylo.sif \
pal2nal "$alignment" "${alignment/_prot.aln/.fasta}" -output fasta > "${alignment/prot.aln/exon_aligned.fasta}"
done
# rename and remove intermediate files
for file in *exon_aligned.fasta; do
mv "$file" "${file/_exon_aligned/}"
rm "${file/exon_aligned.fasta/prot}"*
done
```

## Distances
Convert the overall alignment to distance for generating a network  
Use the GENPOFAD distance model (shouldn't make much difference)
```bash
# from where the full alignment file is
Rscript ~/scripts/align_to_distance.R -p g aligned.fasta
```
Use the resulting `dist_out.nex` to make a network with NeighborNet in SplitsTree4  

As another approach convert the individual gene alignments to distances for generating a network  
Use the GENPOFAD distance model
```bash
# from the "CDS_alignment" folder
Rscript ~/scripts/align_to_distance.R -p g *.fasta
```
Use the resulting `dist_out.nex` to make a network with NeighborNet in SplitsTree4  

## Phylogenetic analysis (partitioned)
To run a maximum likelihood phylogenetic analysis in IQ-TREE, partition by gene and codon position
```bash
# from the "CDS_alignment" folder with 79 genes
python3 ~/scripts/combine_alignments.py -f single -p CDS *.fasta
mkdir mltree && cd mltree
mv ../combine* .
sed -i 's/=/= combine_out.fasta:/g' combine_partitions.nex
singularity exec -H "$(pwd)" ~/singularity-containers/phylo.sif iqtree -T 2 \
--ufboot 1000 --sampling GENESITE -p combine_partitions.nex -m MFP --merge \
--prefix plastome --runs 10
```

Plot the tree if desired  
(use FOR and HERB as outgroups for easier display, as they are on the longest branch)
```bash
grep -E "FOR|HERB" combine_out.fasta | cut -f 2 -d ">" | cut -f 1 -d " " | sort | uniq > outgroup.txt
Rscript ~/scripts/plot_trees.R -b 75 -o outgroup.txt plastome.treefile
mv trees.pdf plastome_tree.pdf
```

# Nuclear ribosomal DNA
Assemble a portion or all of the nuclear ribosomal DNA repeat unit for each sample  
Delimit 18S, ITS1, 5.8S, ITS2 and 26S, and potentially the ETS regions  
Use read mapping to phase differences detected in the ITS and ETS regions  

## Assembly
Given the complexity in the nrDNA repeat regions, it is difficult to assemble them completely with short reads  
Approach:  
1) detect contigs in the SPAdes assemblies that match the nuclear rDNA  
2) pull out longer contigs, align, and generate a consensus as a mapping reference  
3) attempt to assemble for each sample from mapped reads using Unicycler  
4) align assemblies and trim to conserved regions  
5) use mapping to phase sequences when possible  

### Reference generation
Create a new folder `Adenanthos/ribosomal`  

Download references to extract full copies of 18S and 26S rDNA from *Arabidopsis thaliana*:  
18S (partial), 25S, 5.8S	Arabidopsis thaliana	X52320.1  
full repeat	Arabidopsis thaliana	X52322.1  

Also download partial copies from *Adenanthos* and *Isopogon*:  
18S (partial)	Adenanthos sericeus	X66774.1  
26S (partial)	Adenanthos sericeus	X66766.1  
18S (partial)	Isopogon latifolius	X66775.1  
26S (partial)	Isopogon latifolius	X66767.1  

Finally, download portions of ITS from *Adenanthos*:  
ITS1,5.8S,ITS2 (partial)	A. venosus	AF508821.1  
ITS1,5.8S,ITS2 (partial)	A. drummondii	KM659705.1  
ITS1,5.8S,ITS2 (partial)	A. barbiger	AF508822.1  
ITS1,5.8S,ITS2 (partial)	A. linearis	KM659701.1  

Put the accession numbers in a text file (`gb_ref_accessions.txt`), then download
```bash
python3 ~/scripts/get_genbank.py gb_ref_accessions.txt
```

First, concatenate all sequences into a single reference for blasting
```bash
cat *.gb > temp
python3 ~/scripts/genbank_to_fasta.py temp
mv temp.fasta gb_refs.fasta && rm temp*
```

Also pull out copies of *Arabidopsis* 18S, 5.8S and 26S (called "25S") for later use
```bash
# list and extract rRNA from the full repeat
python3 ~/scripts/genbank_parse.py -l rRNA Arabidopsis_thaliana_X52322.gb > temp
python3 ~/scripts/genbank_parse.py -f temp Arabidopsis_thaliana_X52322.gb
# split the extract into 18S and 5.8S
python3 ~/scripts/fasta_splitting.py Arabidopsis_nucl_extract.fasta
mv 18s*.fasta Arabidopsis_18S.fasta
mv 5.8s*.fasta Arabidopsis_58S.fasta
# list and extract rRNA from the portion
python3 ~/scripts/genbank_parse.py -l rRNA Arabidopsis_thaliana_X52320.gb > temp
python3 ~/scripts/genbank_parse.py -f temp Arabidopsis_thaliana_X52320.gb
# split the extract into 5.8S and 25S
python3 ~/scripts/fasta_splitting.py Arabidopsis_nucl_extract.fasta
mv 25s*.fasta Arabidopsis_26S.fasta
rm 5.8s*.fasta Arabidopsis_nucl_extract.fasta temp
```

Upload `gb_refs.fasta` to the same folder on Setonix, then link across the assembly files
```bash
for sample in $(cat ../samples.txt); do
ln -s /scratch/pawsey0220/anderson/Adenanthos/assembly/"$sample"/scaffolds.fasta ./"$sample"_scaffolds.fasta
paste <(echo $(pwd)/${sample}_scaffolds.fasta) <(echo $(pwd)/gb_refs.fasta) >> blast_file.txt
done
```

Run BLAST searches against the scaffolds to detect promising contigs
```bash
cp /software/projects/pawsey0220/anderson/scripts/blast_array.sbatch .
# modify the array script for the number of queries (76, so specify 0-75)
sbatch blast_array.sbatch blast_file.txt
```

Summarise the results
```bash
for sample in $(cat ../samples.txt); do
for line in $(tail -n +2 "$sample"_blast_out.tab | cut -f 1 | sort | uniq); do
paste <(echo $sample) <(echo $line) >> summary.txt
done
done
```

Most samples (46) had multiple contigs hit  
Contig lengths ranged from 683 bp to 16,389 bp  
Pull out the hit scaffolds that are > 6 kbp long and have > 200 coverage
```bash
# generate a list of contigs (43; 2 samples have 2)
awk -F'\t|_' -v OFS='\t' '$5 > 6000 && $7 > 200 {print $1, $2 "_" $3 "_"}' summary.txt > grab_contigs.txt

# use a script to pull out fasta entries; run it in an interactive job
salloc -p work -n 1 -N 1 -c 2 --mem=4G -A pawsey0220 --time=01:00:00
module load python/3.11.6
fasta_script=/software/projects/pawsey0220/anderson/scripts/fasta_extract.py
index=1

while IFS=$'\t' read -r -a myArray; do
sample="${myArray[0]}"
search_string="${myArray[1]}"
python3 "$fasta_script" -s "$search_string" "$sample"_scaffolds.fasta
cat extract.fasta > "$sample"_"$index"_extract.fasta
index="$((index + 1))"
done < grab_contigs.txt

rm extract.fasta
exit
```

Use a reference to re-orient the extracts to the same direction (forward 18S)  
This will also trim them to 500 bp either side of the 18S-ITS-26S piece using the references
```bash
for file in *_extract.fasta; do
bash ~/scripts/reform_ribo.sh "$file" Arabidopsis_18S.fasta Arabidopsis_26S.fasta
if [ -f "new_ribo.fasta" ]; then
mv new_ribo.fasta "${file/_extract/_reform}"
fi
done
```
This resulted in 41 reformed contigs; two contigs were missing hits to 18S or 26S (probably incomplete or off target)  

Align them with MAFFT for inspection in JalView (or other alignment viewer)
```bash
cat *reform.fasta > input.fasta
mafft --auto input.fasta > aligned.fasta
```

There are apparent assembly issues for 11 of the contigs (e.g. inversions or off-target connections)  
Put the NODE names in a `drop_contigs.txt` file, then remove the corresponding fasta files  
```bash
for file in $(grep -f drop_contigs.txt *reform.fasta | cut -f 1 -d ":"); do
rm $file
done
```

Align again with MAFFT, remove sites with less than 50% data, then generate a consensus
```bash
cat *reform.fasta > input.fasta
mafft --auto input.fasta > aligned.fasta
python3 ~/scripts/clean_alignment.py -p 50 aligned.fasta
python3 ~/scripts/consensus.py aligned_clean.fasta
```

Rename the consensus `ribo_reference.fasta` and upload to Setonix for read mapping and assembly  

### Unicycler
Run Unicycler assemblies for all samples, using the `ribo_reference.fasta` as a target for read mapping
```bash
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz .
# remove the full sized HERB files (use the downsampled set of available reads)  
rm HERB{1,2}_R{1,2}.fastq.gz
# create the file with sample names and reference path (one file)
for sample in $(cat ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/ribo_reference.fasta) >> unicycler_file.txt
done
# copy over the scripts
cp /software/projects/pawsey0220/anderson/scripts/unicycler.sh .
# modify the unicycler script to set minid=0.8 and pairlen=800 for lower strictness in mapping
cp /software/projects/pawsey0220/anderson/scripts/unicycler_array.sbatch .
# modify the array script to match how many samples to assemble (76, so 0-75)
sbatch unicycler_array.sbatch unicycler_file.txt
```
The jobs took about 3–30 minutes  

Download the `assembly.gfa` files for each sample  
```bash
for sample in $(cat ../samples.txt); do
cp $sample/unicycler/assembly.gfa ./"$sample"_assembly.gfa
done
# then download the gfa files
```

Use Bandage to open each gfa file, draw the assembly, then blast the two Arabidopsis references (18S and 26S) against it  
In some cases, exporting the assembly will require choosing alternative paths  
Extract the higher depth path for the main copy, as `{sample}_assembly.fasta`  

For samples missing hits to 18S or 26S (too short), these can still be aligned,  
but they may need to be mapped to a different reference for re-assembly/phasing  

In some assemblies, there is a secondary sequence at lower depth with a different hit profile  
If both can be recovered, make a second output with `_b` to check alignment  

If assemblies are fragmented (two or more contigs), also create a second output `_b` for alignment  

Run a script to re-orient the pieces prior to alignment
```bash
for assembly in $(ls *assembly.fasta); do
echo "$assembly"
bash ~/scripts/reform_ribo.sh $assembly Arabidopsis_18S.fasta Arabidopsis_26S.fasta
if [ -f new_ribo.fasta ]; then
sample=$(echo ${assembly/_assembly\.fasta/})
mv new_ribo.fasta "$sample"_reform_assembly.fasta
fi
echo
done
```

Align with MAFFT for inspection  
(also change the headers to sample name rather than contig)
```bash
for assembly in $(ls *reform_assembly.fasta); do
sample=$(echo ${assembly/_reform_assembly\.fasta/})
sed -i "s/^>.*$/>$sample/g" $assembly
done
cat *_reform_assembly.fasta > input.fasta
mafft --auto --thread 8 input.fasta > aligned.fasta
```

The secondary sequences recovered in some samples are clearly divergent  
Running an NCBI BLAST on one of the assemblies hits 93% identity to fungi (Ascomycocta: Pezizomycotina) rDNA  
As this is contamination, the secondary sequences can clearly be ignored  
In two samples, the fungal portion had higher depth reported:  
FOR1-09, FOR1-11  
and in another two, the fungal portion was the only sequence assembled:  
FOR1-04, FOR2-02  

Given the presence of low to high depth fungal contaminants, run a second assembly with stricter mapping (minid = 0.95)  
Use a single representative (non-fungal) assembly from each taxon to map the samples of that taxon to:  
CUN:	CUN1-05  
DOB:	DOB1-01  
FOR:	FOR3-03  
GLA:	GLA2-03  
ILE:	ILE1-02  
SER:	SER2-03  
Map the HERB samples to the FOR reference (there are few differences across the repeat for all taxa)  

Upload the reference assemblies, remove the existing sample folders from the first run,  
and generate a new `unicycler_file.txt` to point the mappings to appropriate references
```bash
rm unicycler_file.txt
for sample in $(grep CUN ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/CUN1-05_reform_assembly.fasta) >> unicycler_file.txt
done
for sample in $(grep DOB ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/DOB1-01_reform_assembly.fasta) >> unicycler_file.txt
done
for sample in $(grep FOR ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/FOR3-03_reform_assembly.fasta) >> unicycler_file.txt
done
for sample in $(grep GLA ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/GLA2-03_reform_assembly.fasta) >> unicycler_file.txt
done
for sample in $(grep HERB ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/FOR3-03_reform_assembly.fasta) >> unicycler_file.txt
done
for sample in $(grep ILE ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/ILE1-02_reform_assembly.fasta) >> unicycler_file.txt
done
for sample in $(grep SER ../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/SER2-03_reform_assembly.fasta) >> unicycler_file.txt
done
```

Change the `unicycler.sh` script to adjust mapping stringency to `minid=0.95`, then run assemblies
```bash
sbatch unicycler_array.sbatch unicycler_file.txt
```
The jobs took about 3.5–37 minutes  

Download the `assembly.gfa` files for each sample  
```bash
for sample in $(cat ../samples.txt); do
cp $sample/unicycler/assembly.gfa ./"$sample"_assembly.gfa
done
# then download the gfa files
```

Use Bandage as before and create `{sample}_assembly.fasta` for each sample  

The following 10 samples only recovered a fragment with 26S:  
CUN2-01, DOB1-02, DOB1-03, FOR1-03, FOR1-14, FOR1-19, FOR3-05, ILE1-04, ILE1-05, ILE2-04  

GLA1-01 and HERB1 were fragmented into two pieces across ITS (save as 18S - `_a` and 26S - `_b`)  

Run a script to re-orient the pieces prior to alignment (remove old assemblies if still present)  
(also correct fasta headers to sample name)
```bash
for assembly in $(ls *assembly.fasta); do
echo "$assembly"
bash ~/scripts/reform_ribo.sh $assembly Arabidopsis_18S.fasta Arabidopsis_26S.fasta
if [ -f new_ribo.fasta ]; then
sample=$(echo ${assembly/_assembly\.fasta/})
sed "s/^>.*$/>$sample/g" new_ribo.fasta > "$sample"_reform_assembly.fasta
rm new_ribo.fasta
fi
echo
done
```

Align with MAFFT for inspection  
```bash
cat *_reform_assembly.fasta > input.fasta
mafft --auto --thread 8 input.fasta > aligned.fasta
```

Based on inspection:  
FOR1-04 - partial fungal in 26S; was detected in Bandage (break in coverage); re-export, excluding that portion  
FOR1-09 - partial fungal in 26S; was detected in Bandage (two options); re-export, choosing *lower* depth option  
FOR2-04 - partial fungal in 26S; not evident in Bandage (slight coverage difference); remove assembly  

Reform the assemblies again, this time only keeping 100 bp on either side of 18S and 26S  
```bash
rm *reform*.fasta FOR2-04_assembly.fasta
for assembly in $(ls *assembly.fasta); do
echo "$assembly"
bash ~/scripts/reform_ribo.sh $assembly Arabidopsis_18S.fasta Arabidopsis_26S.fasta 100
if [ -f new_ribo.fasta ]; then
sample=$(echo ${assembly/_assembly\.fasta/})
sed "s/^>.*$/>$sample/g" new_ribo.fasta > "$sample"_reform_assembly.fasta
rm new_ribo.fasta
fi
echo
done
```

Align with MAFFT for inspection  
```bash
cat *_reform_assembly.fasta > input.fasta
mafft --auto --thread 8 input.fasta > aligned.fasta
```

Based on inspection:  
CUN2-05 has some misassembled bases near the start of 18S  
Except for ITS (missing in HERB1), HERB1 and HERB2 are identical, so the HERB2 assembly could be used for mapping HERB1  

There is some ambiguity for sequences after the end of 26S, so it may be worthwhile trimming references closer than 100 bp  

### Mapping
To inspect assemblies and potentially phase differences, run read mapping to representative sequences for each taxon  

CUN - all samples (where assembled) are identical across the 18S–26S piece; use CUN1-05 for mapping  
DOB - there are a few single bp differences in DOB1-04, but otherwise uniform; use DOB1-05 for mapping  
FOR - there is some ITS variation (substitutions) across the samples, but otherwise mostly uniform; use FOR1-02 for mapping  
GLA - there are a few single bp differences across the samples, but otherwise uniform; use GLA1-05 for mapping  
HERB - HERB2 matches FOR1-02, so both samples can also map against FOR1-02  
ILE - there are a few single bp differences and a single bp indel near ITS, but otherwise uniform; use ILE1-02 for mapping  
SER - SER1-05 has a few single bp difference, but the rest are uniform; use SER1-01 for mapping  

Reform the reference assemblies to keep only 40 bp beyond the ends of 18S and 26S  
Also remove any flanking Ns if present (if there were internal Ns, this step is incorrect!)  
Put the sample names in `ref_samples.txt`
```bash
for sample in $(cat ref_samples.txt); do
echo "$sample"
bash ~/scripts/reform_ribo.sh "$sample"_reform_assembly.fasta Arabidopsis_18S.fasta Arabidopsis_26S.fasta 40
if [ -f new_ribo.fasta ]; then
python3 ~/scripts/remove_ns.py new_ribo.fasta
sed "s/^>.*$/>$sample/g" mod_new_ribo.fasta > "$sample"_reform_trim_assembly.fasta
rm *new_ribo.fasta
fi
echo
done
```

Upload the reference assemblies to a new folder on Setonix `Adenanthos/ribosomal/map`  
In the new folder, create links to the QC reads again, removing the non-downsampled versions of the HERB samples
```bash
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz .
rm HERB{1,2}_R{1,2}.fastq.gz
```

Make a text file with each sample name followed by the reference fasta file path to map to
```bash
for sample in $(grep CUN ../../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/CUN1-05_reform_trim_assembly.fasta) >> mapping_file.txt
done
for sample in $(grep DOB ../../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/DOB1-05_reform_trim_assembly.fasta) >> mapping_file.txt
done
for sample in $(grep FOR ../../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/FOR1-02_reform_trim_assembly.fasta) >> mapping_file.txt
done
for sample in $(grep GLA ../../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/GLA1-05_reform_trim_assembly.fasta) >> mapping_file.txt
done
for sample in $(grep HERB ../../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/FOR1-02_reform_trim_assembly.fasta) >> mapping_file.txt
done
for sample in $(grep ILE ../../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/ILE1-02_reform_trim_assembly.fasta) >> mapping_file.txt
done
for sample in $(grep SER ../../samples.txt); do
paste <(echo ${sample}) <(echo "$(pwd)"/SER1-01_reform_trim_assembly.fasta) >> mapping_file.txt
done
```

Copy over the scripts, adjust for the correct number of samples in the array script, then launch mapping
```bash
cp /software/projects/pawsey0220/anderson/scripts/mapping.sh .
# modify the mapping.sh script to allow mapping of unpaired mates (pairedonly="f")
cp /software/projects/pawsey0220/anderson/scripts/mapping_array.sbatch .
# modify the array script to match how many samples to map (76, so 0-75)
sbatch mapping_array.sbatch mapping_file.txt
```
Took about 30 sec to 3 min per sample  

Run a summary of the mapping using Samtools, and generate a consensus for each sample
```bash
cp /software/projects/pawsey0220/anderson/scripts/samtools_summary.sh .
cp /software/projects/pawsey0220/anderson/scripts/samtools_array.sbatch .
# modify the array script to match how many samples to map (76, so 0-75)
# modify the samtools_summary.sh to use min_depth=20, het_fraction=0.3, call_fraction=0.6
sbatch samtools_array.sbatch -f ../../samples.txt
```
Took about 1 to 39 sec per sample  

Concatenate the summary files together for input to a spreadsheet
```bash
index=1
for sample in $(cat ../../samples.txt); do
if [ $index == 1 ]; then
	head -n 1 "$sample"/"$sample"_samtools_summary.txt > samtools_overall.txt
fi
tail -n 1 "$sample"/"$sample"_samtools_summary.txt >> samtools_overall.txt
index=$((index + 1))
done
```

Align consensus sequences to inspect  
Change the names to match the samples (default is to use the reference name)  
```bash
for sample in $(cat ../../samples.txt); do
sed -i "s/^>.*$/>$sample/g" "$sample"_consensus.fasta
done
```
Also concatenate the *Arabidopsis* fasta references for alignment to help delimit regions  
```bash
cat ../Arabidopsis*.fasta *consensus.fasta > input.fasta
mafft --auto --thread 8 input.fasta > aligned_con.fasta
```

Download the mapping files (*.bam and *.bai and the reference fastas) to visually inspect sites of suspected hets or problems  
Use the Integrated Genomics Viewer or similar program for viewing BAM files and mapping references  

A number of FOR samples had clearly divergent sequences in 18S and 26S that reflect the fungal contamination  
For ease of assessing true heterozygosity, re-map for the FOR samples using more stringent settings (minid = 0.97)  
```bash
grep "FOR" mapping_file.txt | grep -v "HERB" > mapping_file2.txt
# modify the mapping.sh script to make it more strict (minid = 0.97, pairedonly = t)
# modify the array script to match how many samples to map (29, so 0-28)
rm -r FOR*/
sbatch mapping_array.sbatch mapping_file2.txt
```

Run a summary of the mapping using Samtools, and generate a consensus for each sample
```bash
grep "FOR" ../../samples.txt > samples2.txt
# modify the array script to match how many samples to map (29, so 0-28)
sbatch samtools_array.sbatch -f samples2.txt
```

Concatenate the summary files together for input to a spreadsheet
```bash
index=1
for sample in $(cat samples2.txt); do
if [ $index == 1 ]; then
head -n 1 "$sample"/"$sample"_samtools_summary.txt > samtools_overall2.txt
fi
tail -n 1 "$sample"/"$sample"_samtools_summary.txt >> samtools_overall2.txt
index=$((index + 1))
done
```

Rename headers and repeat alignment of new consensus sequences  

Most samples can be phased in either ITS1 and/or ITS2 for minor differences, while still retaining hets in other locations  

For HERB1 and HERB2 (by far the most heterozygous positions), phasing is difficult between ITS1 and ITS2  
The signal in ITS1 is weaker (lower depth), though still evident when examining mapped reads in IGV  
For clarity, focus on phasing the portion of ITS2  
HERB1  
orig:	2305 YTMR..	2362 S..	2434 WWYWSYRAGAKAY  
copy1: 	2305 TTCA..	2362 C..	2434 TATACCAAGAGAT  
copy2: 	2305 CTAG..	2362 G..	2434 ATCTGTGAGATAC  
HERB2  
orig:	2304 YTMR..	2361 S..	2433 WWYWSYRAGAKAY  
copy1:	2304 TTCA..	2361 C..	2433 TATACCAAGAGAT  
copy2:	2304 CTAG..	2361 G..	2433 ATCTGTGAGATAC  

For generating copies for alignment and phylogenetic analysis, use the phased ITS2 (most changes) and  
two positions in ITS1 (so four copies of each HERB sample) that were also variable and phased in the rest of the FOR populations  
That means that other positions in HERB1 and HERB2 remain ambiguities (even though phaseable in places)  
Those positions can be illustrated in a figure, but won't be able to be directly compared in phylogenetic analysis (unphased)  

Generate four text files (`positions1.txt` ..) to create phased copies of the ribosomal region where possible  
Each file should have four columns (no header):  
fasta_id	position	original_bp	new_bp  

Run the base pair altering script on the consensus sequences from the mapping to generate the ribosomal copies  
```bash
for index in 1 2 3 4; do
for sample in $(cut -f 1 positions"$index".txt | sort | uniq); do
python3 ~/scripts/bp_alter.py -f "$sample"_consensus.fasta -p positions"$index".txt
sed "s/${sample}/${sample}_${index}/g" output_alter.fasta > "$sample"_ribo"$index".fasta
rm output_alter.fasta
done
done
```

For those samples not generating phased copies, simply make a copy of the consensus for inclusion  
```bash
ls *ribo*.fasta | cut -f 1 -d _ | sort | uniq > temp
ls *consensus.fasta | cut -f 1 -d _ | sort | uniq > temp2
for sample in $(grep -vf temp temp2); do
cp "$sample"_consensus.fasta "$sample"_ribo.fasta
done
rm temp*
```

Keep the phased copies as well as the consensus sequences for later reference and use in analyses  
It is clear HERB1 and HERB2 are hybrids between FOR and CUN  

## Annotation
To annotate the rRNA portions of nrDNA, Chloe includes a "nr" option
```bash
singularity run -H "$(pwd)" ~/singularity-containers/chloe.sif annotate \
--no-transform -r nr --gbk --embl --no-gff *consensus.fasta
```

For submitting to ENA, there is an online template that could be filled out per sequence  
The suggested format of the flat file indicates a deviation from the normal Chloe output (e.g. no "gene" features per rRNA piece)  
I'm not sure if that means uploading with a "gene" feature will result in rejection  
The ITS pieces are also not annotated by Chloe, although this is suggested for submission  
Note: by default, Chloe converts ambiguous DNA bases to "A", but Chris Jackson created a pull request to correct this (not merged yet)  

Add organism info to the flat files for submission/parsing
```bash
ls *consensus.fasta | cut -f1 -d "_" > samples.txt
for sample in $(grep CUN samples.txt); do
sed -i '1 a SOURCE      nrDNA Adenanthos cuneatus\n  ORGANISM  Adenanthos cuneatus' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos cuneatus"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep DOB samples.txt); do
sed -i '1 a SOURCE      nrDNA Adenanthos dobsonii\n  ORGANISM  Adenanthos dobsonii' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos dobsonii"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep FOR samples.txt); do
sed -i '1 a SOURCE      nrDNA Adenanthos forrestii\n  ORGANISM  Adenanthos forrestii' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos forrestii"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep GLA samples.txt); do
sed -i '1 a SOURCE      nrDNA Adenanthos glabrescens\n  ORGANISM  Adenanthos glabrescens subsp. exasperatus' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos glabrescens subsp. exasperatus"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep HERB samples.txt); do
sed -i '1 a SOURCE      nrDNA Adenanthos eyrei\n  ORGANISM  Adenanthos eyrei' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos eyrei"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep ILE samples.txt); do
sed -i '1 a SOURCE      nrDNA Adenanthos ileticos\n  ORGANISM  Adenanthos ileticos' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos ileticos"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
for sample in $(grep SER samples.txt); do
sed -i '1 a SOURCE      nrDNA Adenanthos sericeus\n  ORGANISM  Adenanthos sericeus' "$sample"_*.gbk
sed -i '4 aFT   source          1..??\nFT                   /organism="Adenanthos sericeus"\nFT                   /mol_type="genomic DNA"' "$sample"_*.embl
done
```

Replace the unknown lengths with the number at the end of the embl files
```bash
for file in *.embl; do
mynum=$(grep Sequence $file | tr -s ' ' | cut -f 3 -d ' ')
sed -i "s/??/${mynum}/" $file
done
```

Replace "circular" with "linear" (note: add spaces for the GenBank header formatting)
```bash
for file in *.gbk; do
sed -i "s/circular/linear  /" $file
done
for file in *.embl; do
sed -i "s/circular/linear/" $file
done
```

Replace Chloe's designation "RRNA" with " rRNA"
```bash
for file in *.gbk; do
sed -i 's/RRNA/ rRNA/g' $file
done
for file in *.embl; do
sed -i 's/RRNA/ rRNA/g' $file
done
```

Replace products (if present), which are in order (18S, 5.8S, 26S) and "missing" given my pull request (or completely absent with base Chloe)
```bash
for file in *.gbk; do
sed -i '0,/missing/{s/missing/18S ribosomal RNA/}' $file
sed -i '0,/missing/{s/missing/5.8S ribosomal RNA/}' $file
sed -i '0,/missing/{s/missing/26S ribosomal RNA/}' $file
done
for file in *.embl; do
sed -i '0,/missing/{s/missing/18S ribosomal RNA/}' $file
sed -i '0,/missing/{s/missing/5.8S ribosomal RNA/}' $file
sed -i '0,/missing/{s/missing/26S ribosomal RNA/}' $file
done
```

Remove extraneous lines from the GenBank and ENA files
```bash
for sample in $(cat samples.txt); do
sed -i '/ID=/d' "$sample"_*.gbk
sed -i '/ID=/d' "$sample"_*.embl
sed -i '/Parent/d' "$sample"_*.gbk
sed -i '/Parent/d' "$sample"_*.embl
sed -i '/Name/d' "$sample"_*.gbk
sed -i '/Name/d' "$sample"_*.embl
sed -i '/locus_tag/d' "$sample"_*.gbk
sed -i '/locus_tag/d' "$sample"_*.embl
done
```

Annotate the ITS pieces as "misc_RNA", using four numbers for the boundaries
```bash
for sample in $(cat samples.txt); do
file="$sample"_consensus.chloe.embl
# grab coordinates for end of 18S, start of 5.8S, end of 5.8S, start of 26S
end18s=$(grep "18S rRNA" -B 1 $file | head -n 1 | cut -f 3 -d ".")
start58s=$(grep "5.8S rRNA" -B 1 $file | head -n 1 | tr -s ' ' | cut -f 3 -d ' ' | cut -f 1 -d ".")
end58s=$(grep "5.8S rRNA" -B 1 $file | head -n 1 | cut -f 3 -d ".")
start26s=$(grep "26S rRNA" -B 1 $file | head -n 1 | tr -s ' ' | cut -f 3 -d ' ' | cut -f 1 -d ".")
# insert the ITS features
sed -i "12 aFT   misc_RNA        $((end18s + 1))..$((start58s - 1))\nFT                   /note=\"internal transcribed spacer 1, ITS1\"" "$file"
sed -i "19 aFT   misc_RNA        $((end58s + 1))..$((start26s - 1))\nFT                   /note=\"internal transcribed spacer 2, ITS2\"" "$file"
file="$sample"_consensus.chloe.gbk
sed -i "9 a\ \ \ \ \ misc_RNA        $((end18s + 1))..$((start58s - 1))\n                     /note=\"internal transcribed spacer 1, ITS1\"" "$file"
sed -i "16 a\ \ \ \ \ misc_RNA        $((end58s + 1))..$((start26s - 1))\n                     /note=\"internal transcribed spacer 2, ITS2\"" "$file"
done
```

## Distances
Align the phased copies (ambiguities still present) and estimate genetic distances  
```bash
# from the "map" folder with the consensus sequences and ribo copies
cat *ribo*.fasta > input.fasta
mafft --thread 8 --auto input.fasta > aligned.fasta
Rscript ~/scripts/align_to_distance.R -p g aligned.fasta
```
The output `dist_out.nex` can be use to generate a network using NeighborNet in SplitsTree4  

## Phylogenetic tree (full alignment)
Use IQ-TREE to run a maximum likelihood phylogenetic analysis with 1000 ultrafast bootstraps (UFB) and search for best model  
```bash
# in a new folder "mltree", after copying "aligned.fasta" into it
singularity exec -H "$(pwd)" ~/singularity-containers/phylo.sif iqtree -T 1 \
--ufboot 1000 -m MFP --prefix ribo -s aligned.fasta --runs 10
```

Plot the tree  
(root on SER as a divergent and longer branch)
```bash
grep "SER" aligned.fasta | cut -f 2 -d ">" | sort | uniq > outgroup.txt
Rscript ~/scripts/plot_trees.R -b 75 -o outgroup.txt ribo.treefile
```

# Nuclear genes
It may be possible with sufficient read depth to recover nuclear genes for comparison  
One way to approach this is to use target capture pipelines with the raw data  
Start by assembling from the samples with the most reads (the herbarium putative hybrids)  

Two capture pipelines are HybPiper and Captus  

## HybPiper for A353
One set of putatively low copy genes are the Angiosperms353 genes (A353)  

Use the HybPiper pipeline v. 2.3.2 (https://github.com/mossmatters/HybPiper) to assemble them  

### Target file
Augment the default target file with more Proteaceae references  

Use the "NewTargets" (https://github.com/chrisjackson-pellicle/NewTargets) approach with their `mega353.fasta` file  
Run the script `filter_mega353.py`, with a config file (`select_file.txt`) with the lines  
```
[Family]
Proteaceae
```
```bash
unzip mega353.fasta.zip
python3 filter_mega353.py mega353.fasta select_file.txt
```
The resulting output was named `targets_Proteaceae.fasta`  

Translate the targets to protein  
```bash
python3 ~/scripts/translate.py targets_Proteaceae.fasta
```
Remove any "*" stop codons at the end of sequences; two sequences had internal stop codons and were removed  
Use the resulting `targets_Proteaceae_prot.fasta` for assembly  
Upload it to the `Adenanthos` folder  

#### Preliminary run for HERB
Run assemblies with HERB1 and HERB2 to use the recovered sequences as a new target file for an assembly of all samples  

In a `hybpiper_HERB` folder, run with all the reads available
```bash
mkdir reads
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/HERB*/*.gz reads/
rm reads/*down*		# get rid of the downsampled
echo "HERB1" > samples.txt
echo "HERB2" >> samples.txt
cp ../targets_Proteaceae_prot.fasta .
cp /software/projects/pawsey0220/anderson/scripts/hybpiper_array.sbatch .
# manually change to 2 samples (0-1)
# also double the memory per job to 128G
sbatch hybpiper_array.sbatch -f targets_Proteaceae_prot.fasta -r reads/ -s samples.txt -t "aa"
```
The runs took about 45 minutes  

After the runs finish, combine results
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/combine_hybpiper.sbatch -d "yes"
```

Assess the output to make a new target file  
```bash
# check for loci with retained paralogs in the "results" folder
grep ">" paralog_seqs/*.fasta | grep ".main" | cut -f 1 -d "." | cut -f 2 -d "/" | sort | uniq > paralog_loci.txt
#  11 loci
# remove loci with excessive Ns (incomplete) and shorter than 300 bp
python3 ~/scripts/seq_stats.py dna_seqs/*.fasta > temp
cat temp | grep -v Total | grep "fasta" | awk '$5 < 20 && $3 > 300' | cut -f 2 -d "/" | cut -f 1 -d "." | sort | uniq > loci.txt
# exclude the paralogs
grep -vf paralog_loci.txt loci.txt > keep_loci.txt
# leaving 291 loci
# generate consensus sequences of the two HERB samples per locus
for locus in $(cat keep_loci.txt); do
mafft --auto dna_seqs/"$locus".fasta > "$locus".aln
python3 ~/scripts/consensus.py "$locus".aln
sed "s/^>consensus/>HERB-${locus}/" consensus.fasta >> HERB_targets_ref.fasta
rm "$locus".aln consensus.fasta
done
# translate to protein residues
python3 ~/scripts/translate.py HERB_targets_ref.fasta
# remove any loci with internal stop codons "*"
# put locus names to remove in a "temp" file after finding the "*" and noting the loci
python3 ~/scripts/remove_fastas.py -f temp HERB_targets_ref*.fasta
for file in mod*; do mv $file ${file/mod_/}; done
# this leaves 284 loci in the target files
```

The resulting file `HERB_targets_ref_prot.fasta` can be used as a target file for assembly of all samples  

### Assembly
Make a new folder `hybpiper` with the targets file in it and add a `samples.txt` file with all sample names  
```bash
mkdir reads
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz reads/
rm reads/*down*		# get rid of the downsampled HERB samples
cp ../samples.txt .
cp ../HERB_targets_ref_prot.fasta .
cp /software/projects/pawsey0220/anderson/scripts/hybpiper_array.sbatch .
# manually change the array script to 76 samples (0-75) and double the memory per job to 128G
sbatch hybpiper_array.sbatch -f HERB_targets_ref_prot.fasta -r reads/ -s samples.txt -t "aa"
```
Most samples took about 8 to 20 minutes  
Two of the larger runs took about 40 to 60 minutes  

After the runs finish, combine results
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/combine_hybpiper.sbatch -d "yes"
```
Download the `results/combined_stats.tsv` for assessment  

HybPiper locus recovery was variable across samples, with 27 samples having < 50 loci at > 50% target length  
To extract hopefully reliable loci for alignment, drop undesirable samples and loci with too much missing data  
Remove samples with fewer than 120 loci recovered at all (equivalent to < 50 loci at > 50% target length)  
Remove loci with fewer than 75% of the samples remaining (49, so threshold of c. 37)
```bash
# from the "results" folder
ls dna_seqs/*.fasta | cut -f 2 -d "/" | cut -f 1 -d "." > loci.txt		# 276 loci
# determine samples with few loci
for sample in $(cat ../../samples.txt); do
paste <(echo $sample) <(grep ">${sample}" dna_seqs/*.fasta | wc -l) >> temp
done
awk '$2 < 120' temp | cut -f 1 > sample_drops.txt		# 27 samples
rm temp
# determine loci with few samples (other than excluded ones)
for locus in $(cat loci.txt); do
paste <(echo $locus) <(grep ">" dna_seqs/${locus}.fasta | grep -vf sample_drops.txt | wc -l) >> temp
done
awk '$2 < 37' temp | cut -f 1 > locus_drops.txt		# 114 loci
rm temp
```
Based on this first quick pass, there should be 49 samples and 162 candidate loci for alignment  

### Exons
First, align candidate exons for assessment of differences  
From the `hybpiper` directory, make a new directory `filtered_exons` and change into it
```bash
mkdir filtered_exons && cd filtered_exons
ln -s ../results/dna_seqs/*.fasta .
for droplocus in $(cat ../results/locus_drops.txt); do
rm "$droplocus".fasta
done
# get rid of extra text in the files and copy them over the links
for file in *.fasta; do
sed -E '/>/s/\s(.+)$//g' "$file" > temp && mv temp $file
done
```

Make a subfolder `aligned` and run translations and alignments of protein residues, then back-substitution to nucleotides  
For each file, remove the low cover samples (if present), translate to protein, substitute stop codons "*" with "X",  
align with MAFFT, then use `pal2nal.pl` (Suyama et al. 2006; http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz)  
to convert the nucleotide files to align with the protein alignment
```bash
mkdir aligned && cd aligned
sbatch /software/projects/pawsey0220/anderson/scripts/translate_align.sbatch -f ../ -d ../../results/sample_drops.txt
```

Manually review the alignments for obvious problems and paralogy  
Because of the way HybPiper assembles loci, paralogy can be apparent in single exons but not in others  
When taxa are split into two major copy types (e.g. half the CUN samples have one or the other),  
it is likely because that exon is paralogous  
Based on the review, 57 loci had evident problems or showed signs of paralogy (put these in a `drop_loci.txt` file)  
Others had problematic sequences which could be removed follow initial tree generation  

Remove the paralogous/bad loci
```bash
for locus in $(cat drop_loci.txt); do
rm $locus.fasta
done
```
This leaves 105 loci to analyse  

Run an initial phylogenetic analysis of the loci, followed by TreeShrink to remove problematic sequences from each locus
```bash
cd ..
mkdir phylo && cd phylo
sbatch /software/projects/pawsey0220/anderson/scripts/align_phylo.sbatch -a "loci" -c "n" -con "n" -f "../aligned" -r "n"
```

Assess branch lengths for the loci using TreeShrink and drop outlier samples from the `filtered_exons/aligned` files
```bash
# start an interactive job from the `filtered_exons/phylo` directory
salloc -p work -n 1 -N 1 -c 2 --mem=4G -A pawsey0220 --time=01:00:00
module load singularity/4.1.0-slurm python/3.11.6

# run assessment
singularity exec -H "$(pwd)" /software/projects/pawsey0220/anderson/singularity-containers/phylo.sif run_treeshrink.py \
-t loci.treefile -m per-gene -q 0.10 -O output_ts -o treeshrink

# iterate through the loci, dropping samples detected by TreeShrink
ls ../aligned/*.fasta | cut -f 3 -d "/" | cut -f 1 -d "." > loci.txt
index=1
for locus in $(cat loci.txt); do
remove_line=$(sed -n "${index}p; $((index + 1))q" treeshrink/output_ts.txt | tr -s '\t' ',' | sed 's/,$//')
if [ ! -z "$remove_line" ]; then
python3 /software/projects/pawsey0220/anderson/scripts/remove_fastas.py "$remove_line" ../aligned/"${locus}.fasta"
else
echo "locus $locus does not have taxa to remove"
fi
((index+=1))
done
mkdir shrunk
cp ../aligned/*.fasta shrunk/
for file in mod_*; do
mv "$file" shrunk/"${file/mod_/}"
done

# end the job
exit
```

Download the `shrunk` folder alignments for another check for any sequences that need to be removed (misassembled portions)  

4932:	CUN1-01,FOR1-17,FOR1-20  
4951:	CUN1-01,FOR1-16,FOR1-20,ILE2-01,ILE2-03,SER1-05,SER2-04
4992:	drop (multiple samples with divergent linked differences: probably some paralogs)  
5264:	CUN1-05,DOB1-05,FOR4-01,GLA2-04  
5296:	CUN2-03,FOR1-15,FOR1-21  
5449:	GLA2-04  
5614:	FOR1-17,ILE2-01  
5816:	DOB1-04,SER2-04  
5821:	CUN2-03,FOR1-14,FOR1-16,FOR1-17,FOR1-20  
5921:	FOR1-14,FOR1-16,FOR1-17,FOR1-20,FOR1-21,ILE2-01,SER1-03,SER1-05,SER2-01,SER2-04  
5926:	SER1-03  
6139:	FOR4-01,ILE1-03  
6295:	FOR1-14,FOR1-15,FOR1-16,FOR1-17,FOR1-20,FOR1-21,FOR3-04,FOR4-01  
6363:	FOR1-15,FOR1-20  
6487:	FOR1-15,FOR1-17,SER2-01  
7361:	FOR1-16,ILE1-02  

```bash
# start an interactive job from the `filtered_exons/phylo/shrunk` directory
salloc -p work -n 1 -N 1 -c 2 --mem=4G -A pawsey0220 --time=01:00:00
module load python/3.11.6
# manually remove sequences using a script, e.g.
python3 /software/projects/pawsey0220/anderson/scripts/remove_fastas.py CUN1-01,FOR1-17,FOR1-20 4932.fasta
# (repeat for the others above)
for file in mod_*; do
mv "$file" "${file/mod_/}"
done
exit
```

Re-run the analysis on the cleaned 104 loci in a new folder `phylo/run2`
```bash
mkdir run2 && cd run2
sbatch /software/projects/pawsey0220/anderson/scripts/align_phylo.sbatch -a "all" -c "n" -con "n" -f "../shrunk" -r "n" -p "y"
```

Run concordance analysis for completeness
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/align_phylo.sbatch -a "none" -c "n" -con "y" -f "in_align" -r "n"
```

Calculate summary statistics on the alignments  
```bash
python3 ~/scripts/seq_stats.py in_align/*.fasta > stats.txt
```

Optionally plot trees interactively with `plot_trees.rmd`; need an `outgroup.txt` file  
Put all the SER sample names in an `outgroup.txt` file  
Convert concordance output and astral output using scripts for reading into the plotting markdown  
```bash
python3 ~/scripts/concord_to_newick.py -t concord_gcf.cf.tree.nex -o concord_newick
python3 ~/scripts/concord_to_newick.py -t concord_scf.cf.tree.nex -o concord_newick
python3 ~/scripts/astral_parse.py -t astral.tre -f p -o astral
python3 ~/scripts/astral_parse.py -t astral.tre -f q -o astral
# collapse branches with ASTRAL polytomy p-values >0.05
nw_ed astral_poly.tre "i & b > 0.05" o > astral_poly_collapsed.tre
```
Plot the trees interactively with `plot_trees.rmd`  

Alternatively, plot simpler trees without concordance with `plot_trees.R`
```bash
Rscript ~/scripts/plot_trees.R -b 75 -o outgroup.txt concat.treefile && mv trees.pdf concat_simple.pdf
# the astral posterior probability is the first of number/number/number on branch labels
# reduce numbers of decimal places
sed -E 's/([0-9]+\.[0-9][0-9][0-9])[0-9]+/\1/g' astral_p.tre > temp
# grab the first element (possibly some scientific notation)
sed -E 's/([0-9]+\.[0-9]+[E-]*[0-9]*)\/[0-9]+\.[0-9]+[E-]*[0-9]*\/[0-9]+\.[0-9]+[E-]*[0-9]*:/\1:/g' temp > temp2
Rscript ~/scripts/plot_trees.R -b 0.8 -o outgroup.txt temp2 && mv trees.pdf astral_simple.pdf
rm temp*
```

### Supercontigs
Use the selected/clean exon alignments to retrieve associated supercontigs (exons + introns) from the same samples  
These may be useful for mapping reads and consensus generation to recover heterozygosity  

Make a new directory `hybpiper/filtered_super`
```bash
mkdir filtered_super && cd filtered_super
ls ../filtered_exons/phylo/run2/in_align/*.fasta | cut -f 6 -d "/" | cut -f 1 -d "." > loci.txt
for locus in $(cat loci.txt); do
ln -s /scratch/pawsey0220/anderson/Adenanthos/hybpiper/results/supercontig_seqs/"$locus".fasta .
done
```

For each locus, create a list of samples to drop (not matching the ones run in the exon analysis), then drop them
```bash
# start an interactive job
salloc -p work -n 1 -N 1 -c 4 -A pawsey0220 --time=01:00:00
module load python/3.11.6
for locus in $(cat loci.txt); do
grep ">" ../filtered_exons/phylo/run2/in_align/"$locus".fasta | cut -f 2 -d ">" > samples_"$locus".txt
grep -vf samples_"$locus".txt ../../samples.txt > temp_drop.txt
python3 /software/projects/pawsey0220/anderson/scripts/remove_fastas.py -f temp_drop.txt "$locus".fasta
mv mod_"$locus".fasta "$locus".fasta
done
rm samples_*.txt temp*
exit
```

Trim fasta descriptions and format the supercontig sequences to increase Ns spacing to 100  
(Note: may sometimes increase existing strings of Ns to longer, but should be OK for mapping)
```bash
for file in *.fasta; do
sed -E '/>/s/\s(.+)$//g' "$file" > temp && mv temp $file
done
salloc -p work -n 1 -N 1 -c 4 -A pawsey0220 --time=01:00:00
module load python/3.11.6
python3 /software/projects/pawsey0220/anderson/scripts/fasta_Ns_sub.py *.fasta
exit
```

For each locus, split the individual fastas out, then move them to their own sample directories
```bash
# start an interactive job
salloc -p work -n 1 -N 1 -c 4 -A pawsey0220 --time=01:00:00
module load python/3.11.6
for locus in $(cat loci.txt); do
python3 /software/projects/pawsey0220/anderson/scripts/fasta_splitting.py "$locus".fasta
for sample in $(grep ">" "$locus".fasta | cut -f 2 -d ">"); do
if [ ! -d "$sample" ]; then
mkdir "$sample"
fi
if [ -f "$sample".fasta ]; then
mv "$sample".fasta "$sample"/"$locus".fasta
else
echo "Missing file for $sample for locus $locus"
fi
done
done
# include the locus number in the fasta header for each locus for each sample
grep ">" *.fasta | cut -f 2 -d ">" | sort | uniq > map_samples.txt
for sample in $(cat map_samples.txt); do
cd $sample
for file in *.fasta; do
locus=${file/\.fasta/}
sed -i "/^>/s/$/-$locus/" $file
done
cd ..
done
exit
```

Link reads, concatenate the loci together per sample, map, then generate a consensus per sample
```bash
mkdir mapping && cd mapping
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz .
rm HERB*_down_*.gz		# remove the downsampled HERB reads

# concatenate loci and make mapping file
for sample in $(cat ../map_samples.txt); do
cat ../"$sample"/*.fasta > "$sample"_map_ref.fa
paste <(echo ${sample}) <(echo ${sample}_map_ref.fa) >> mapping_file.txt
done

# copy over scripts and run mapping
cp /software/projects/pawsey0220/anderson/scripts/mapping.sh .
# modify the mapping script to make it more strict (minid=0.95)
cp /software/projects/pawsey0220/anderson/scripts/mapping_array.sbatch .
# modify the array script to match how many samples to map (49, so 0-48)
sbatch mapping_array.sbatch mapping_file.txt

# after it finishes, copy over scripts for samtools and run consensus generation
cp /software/projects/pawsey0220/anderson/scripts/samtools_summary.sh .
# modify the samtools_summary.sh to use min_depth=6, het_fraction=0.3, call_fraction=0.6
cp /software/projects/pawsey0220/anderson/scripts/samtools_array.sbatch .
# modify the array script to match how many samples to map (49, so 0-48)
sbatch samtools_array.sbatch -f ../map_samples.txt
```

Concatenate the summary files together for input to a spreadsheet
```bash
index=1
for sample in $(cat ../map_samples.txt); do
if [ $index == 1 ]; then
head -n 1 "$sample"/"$sample"_samtools_summary.txt > samtools_overall.txt
fi
tail -n 1 "$sample"/"$sample"_samtools_summary.txt >> samtools_overall.txt
index=$((index + 1))
done
# NOTE: the positions of hets and Ns will be extensive, so could drop those columns for limiting file size
# cut -f 8,10 --complement samtools_overall.txt > samtools_overall_smaller.txt
```

Split out the fastas from the consensus sequences and move to each sample directory
```bash
# start an interactive job
salloc -p work -n 1 -N 1 -c 4 -A pawsey0220 --time=01:00:00
module load python/3.11.6

# grab the consensus sequences and move them
for sample in $(cat ../map_samples.txt); do
python3 /software/projects/pawsey0220/anderson/scripts/fasta_splitting.py "$sample"/"$sample"_consensus.fasta
for outfasta in "$sample"-*.fasta; do
mv "$outfasta" ../"$sample"/"${outfasta/${sample}-/consens_}"
done
done

# change back the "filtered_super" folder and concatenate across samples per locus, removing locus info from headers
cd ..
mkdir aligning
for locus in $(cat loci.txt); do
for sample in $(cat map_samples.txt); do
if [ -f "$sample"/consens_"$locus".fasta ]; then
sed "/^>/s/-$locus//" "$sample"/consens_"$locus".fasta >> aligning/"$locus".fasta
fi
done
done
exit
```

Run alignments of the loci with MAFFT
```bash
cp /software/projects/pawsey0220/anderson/scripts/align_phylo.sbatch .
# modify the script for the cleaning step to remove positions with >80% gaps and keep samples regardless of missing
sbatch align_phylo.sbatch -a "none" -c "y" -f "aligning" -r "y"
```

Generate a distance matrix from the alignments for use with SplitsTree4
```bash
Rscript ~/scripts/align_to_distance.R -p "g" in_align/*.fasta
```
The resulting distance matrix is not well resolved nor particularly informative  
It is also not substantially more resolved or different than a similar distance network based on the exons only  
HERB1 and HERB2 are more closely associated and "nearer" to SER and CUN, but not clearly  

The consensus sequences and ambiguities may be poorly recovered or uninformative at heterozygous positions(?)  

Heterozygosity is nonetheless clearly different (outliers) for the HERB samples:  
HERB samples: 1.13% heterozygosity (1.05–1.21)  
other samples: 0.53% heterozygosity (0.3–0.7)  

## Captus
As an alternative to the HybPiper pipeline, Captus performs an assembly of the cleaned sequences first before pulling out loci  

Make a directory `captus` in the top `Adenanthos` directory and change into it  
Copy over the same target file as used for the HybPiper run (based on initial run for HERB samples)
```bash
mkdir reads
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz reads/
rm reads/*down*		# get rid of the downsampled HERB samples
cp ../hybpiper/HERB_targets_ref_prot.fasta .
```

Launch the assembly step (step 2) of Captus  
Modify the script to use 32 cores per sample (4 samples at a time)  
Because the samples have many reads, assemblies take a while, so split the process up into sets of 8 samples
```bash
cp /software/projects/pawsey0220/anderson/scripts/captus.sbatch .
# manually modify to divide cores by 32

# create the sets of samples
ls reads/*.gz > temp
split temp set_ -l 16
for dataset in set_*; do
mkdir -p fold_"$dataset"/reads
for file in $(cat $dataset); do
ln -s "$(pwd)"/$file fold_"$dataset"/reads/
done
done

# run the jobs
for folder in fold_*; do
cd $folder
sbatch ../captus.sbatch -s 2 -r reads
cd ..
done
rm set_*
rm temp
```
Assemblies took about 1 to 3 hours per sample, so jobs of 8 samples took about 6 to 8 hours  
There is no evident improvement in processing time using MEGAHIT compared to SPAdes  

After the jobs finish, move all the files into the main directory in a new folder `02_assemblies`
```bash
mkdir 02_assemblies
mv fold_*/*.log 02_assemblies/
# start an interactive job to copy files
salloc -p copy -n 1 -N 1 -c 2 --mem=8G -A pawsey0220 --time=02:00:00
index=1
for folder in fold_*; do
for file in "$folder"/02_assemblies/*.tsv; do
name=$(basename $file)
mv "$file" 02_assemblies/"${name/stats/stats${index}}"
done
mv "$folder"/02_assemblies/captus-assemble_report.html 02_assemblies/captus-assemble_report${index}.html
rsync -vr --progress "$folder"/02_assemblies/* 02_assemblies/
index=$((index + 1))
done
# remove the folders
find -P fold_* -type f -print0 -o -type l -print0 | xargs -0 munlink
find -P fold_* -type d -empty -delete
# end the job
exit
```
Download the HTML reports to assess assembly success and metrics  

Launch the extraction step (step 3)  
An earlier run to look for additional markers (mode "e") failed to complete in the time limit so don't run
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/captus.sbatch -s 3 -t HERB_targets_ref_prot.fasta -m n
```
Download the HTML report to assess recovery  
A number of loci have somewhat lower average weighted score (potentially bad fits?):  
5703,5464,5354,6405,5554,6563,5421,6540,7021  
Leave them in  

Download the `captus-extract_stats.tsv` file to summarise the recovery
```bash
# fields of interest: 1 - sample, 3 - locus, 10 - pct_identity
tail -n +4 captus-extract_stats.tsv | cut -f 1,3,10 > temp
cut -f 1 temp | sort | uniq > extract_samples.txt		# 76
for sample in $(cat extract_samples.txt); do
grep "$sample" temp > temp2
paste <(echo $sample) <(cut -f 2 temp2 | sort | uniq | wc -l) \
<(cut -f 2 temp2 | sort | uniq -cd | wc -l) \
<(awk '{cnt[$2]++; val[$2] = val[$2] + $3} END {for (key in val) {sum += val[key] / cnt[key]; num++}} END {print sum / num}' temp2) \
>> summary.txt
done
```
There is substantially more reported recovery across samples than with HybPiper  
A number of samples had lower recovery and should probably be dropped (< 250 loci):  
FOR1-09,FOR1-04,FOR1-03,FOR1-11,FOR2-04,FOR1-19,FOR3-03  
Simply remove their `*captus-ext` folders in `03_extractions`  

Run the alignment and trimming (step 4)
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/captus.sbatch -s 4
```
Download the HTML report to assess alignment completeness  
A few samples have more evident gaps and lower completeness:  
FOR1-07,FOR3-01,FOR1-02,FOR1-06,FOR4-01  
Leave them in  

### Exons
Run a phylogenetic analysis of the protein-aware alignment of nucleotides for the exons  
From the `Adenanthos/captus` folder, create new folders `phylo/temp_input` to prepare to run the analysis
```bash
mkdir -p phylo/temp_input && cd phylo
ln -s /scratch/pawsey0220/anderson/Adenanthos/captus/04_alignments/03_trimmed/05_naive/01_coding_NUC/02_NT/*.fna temp_input/
for file in temp_input/*.fna; do
mv $file ${file/\.fna/\.fasta}
done
```

Run an initial analysis of the loci to enable assessment with TreeShrink
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/align_phylo.sbatch -a "loci" -c "n" -con "n" -f "temp_input" -r "n"
```

Run assessment with TreeShrink to detect potential misassemblies or problems
```bash
# start an interactive job from the `filtered_exons/phylo` directory
salloc -p work -n 1 -N 1 -c 2 --mem=4G -A pawsey0220 --time=01:00:00
module load singularity/4.1.0-slurm python/3.11.6

# run assessment
singularity exec -H "$(pwd)" /software/projects/pawsey0220/anderson/singularity-containers/phylo.sif run_treeshrink.py \
-t loci.treefile -m per-gene -q 0.10 -O output_ts -o treeshrink

# iterate through the loci, dropping samples detected by TreeShrink
ls in_align/*.fasta | cut -f 2 -d "/" | cut -f 1 -d "." > loci.txt
index=1
for locus in $(cat loci.txt); do
remove_line=$(sed -n "${index}p; $((index + 1))q" treeshrink/output_ts.txt | tr -s '\t' ',' | sed 's/,$//')
if [ ! -z "$remove_line" ]; then
python3 /software/projects/pawsey0220/anderson/scripts/remove_fastas.py "$remove_line" in_align/"${locus}.fasta"
else
echo "locus $locus does not have taxa to remove"
fi
((index+=1))
done
mkdir shrunk
cp in_align/*.fasta shrunk/
for file in mod_*; do
mv "$file" shrunk/"${file/mod_/}"
done

# end the job
exit
```

Run a second analysis using the adjusted alignment files
```bash
mkdir run2 && cd run2
sbatch --ntasks=2 --cpus-per-task=32 /software/projects/pawsey0220/anderson/scripts/align_phylo.sbatch -a "all" -c "n" -con "n" -f "../shrunk" -r "n" -p "y"
```

Run concordance analysis for completeness
```bash
sbatch /software/projects/pawsey0220/anderson/scripts/align_phylo.sbatch -a "none" -c "n" -con "y" -f "in_align" -r "n"
```

Calculate summary statistics on the alignments  
```bash
python3 ~/scripts/seq_stats.py in_align/*.fasta > stats.txt
```

Optionally plot trees interactively with `plot_trees.rmd`; need an `outgroup.txt` file  
Put all the SER sample names in an `outgroup.txt` file  
Convert concordance output and astral output using scripts for reading into the plotting markdown  
```bash
python3 ~/scripts/concord_to_newick.py -t concord_gcf.cf.tree.nex -o concord_newick
python3 ~/scripts/concord_to_newick.py -t concord_scf.cf.tree.nex -o concord_newick
python3 ~/scripts/astral_parse.py -t astral.tre -f p -o astral
python3 ~/scripts/astral_parse.py -t astral.tre -f q -o astral
# collapse branches with ASTRAL polytomy p-values >0.05
nw_ed astral_poly.tre "i & b > 0.05" o > astral_poly_collapsed.tre
```
Plot the trees interactively with `plot_trees.rmd`  

Alternatively, plot simpler trees without concordance with `plot_trees.R`
```bash
Rscript ~/scripts/plot_trees.R -b 75 -o outgroup.txt concat.treefile && mv trees.pdf concat_simple.pdf
# the astral posterior probability is the first of number/number/number on branch labels
# reduce numbers of decimal places
sed -E 's/([0-9]+\.[0-9][0-9][0-9])[0-9]+/\1/g' astral_p.tre > temp
# grab the first element (possibly some scientific notation)
sed -E 's/([0-9]+\.[0-9]+[E-]*[0-9]*)\/[0-9]+\.[0-9]+[E-]*[0-9]*\/[0-9]+\.[0-9]+[E-]*[0-9]*:/\1:/g' temp > temp2
Rscript ~/scripts/plot_trees.R -b 0.8 -o outgroup.txt temp2 && mv trees.pdf astral_simple.pdf
rm temp*
```

### Genes (exons + introns)
Run read mapping to the full gene sequences recovered to attempt to incorporate heterozygosity  

Make a new consensus directory in the top `captus` directory and change to it  
Link the untrimmed alignments and clean up headers
```bash
mkdir consensus && cd consensus
ln -s /scratch/pawsey0220/anderson/Adenanthos/captus/04_alignments/02_untrimmed/05_naive/01_coding_NUC/03_genes/*.fna .
ls *.fna | cut -f 1 -d "." > loci.txt
for locus in $(cat loci.txt); do
sed -E '/>/s/\s(.+)$//' "$locus".fna > "$locus".fasta
done
rm *.fna
```

For each locus, split the individual fastas out, remove gaps from sequences, then move them to their own sample directories
```bash
# start an interactive job
salloc -p work -n 1 -N 1 -c 4 -A pawsey0220 --time=01:00:00
module load python/3.11.6
for locus in $(cat loci.txt); do
python3 /software/projects/pawsey0220/anderson/scripts/remove_ns.py -r gaps "$locus".fasta
python3 /software/projects/pawsey0220/anderson/scripts/fasta_splitting.py mod_"$locus".fasta
for sample in $(grep ">" "$locus".fasta | cut -f 2 -d ">"); do
if [ ! -d "$sample" ]; then
mkdir "$sample"
fi
if [ -f "$sample".fasta ]; then
mv "$sample".fasta "$sample"/"$locus".fasta
else
echo "Missing file for $sample for locus $locus"
fi
done
rm mod_*.fasta
done
# include the locus number in the fasta header for each locus for each sample
grep ">" *.fasta | cut -f 2 -d ">" | sort | uniq > map_samples.txt
for sample in $(cat map_samples.txt); do
cd $sample
for file in *.fasta; do
locus=${file/\.fasta/}
sed -i "/^>/s/$/-$locus/" $file
done
cd ..
done
exit
```

Make a directory to do the mapping, link reads, concatenate the loci together per sample, map, then generate a consensus per sample
```bash
mkdir mapping && cd mapping
ln -s /scratch/pawsey0220/anderson/Adenanthos/qc/*/*.gz .
rm HERB*_down_*.gz		# remove the downsampled HERB reads (or keep them for comparable mapping, but change labels appropriately)

# concatenate loci and make mapping file
for sample in $(cat ../map_samples.txt); do
cat ../"$sample"/*.fasta > "$sample"_map_ref.fa
paste <(echo ${sample}) <(echo ${sample}_map_ref.fa) >> mapping_file.txt
done

# copy over scripts and run mapping
cp /software/projects/pawsey0220/anderson/scripts/mapping.sh .
# modify the mapping script to make it more strict (minid=0.95)
cp /software/projects/pawsey0220/anderson/scripts/mapping_array.sbatch .
# modify the array script to match how many samples to map (69, so 0-68)
sbatch mapping_array.sbatch mapping_file.txt

# after it finishes, copy over scripts for samtools and run consensus generation
cp /software/projects/pawsey0220/anderson/scripts/samtools_summary.sh .
# modify the samtools_summary.sh to use min_depth=6, het_fraction=0.3, call_fraction=0.6
cp /software/projects/pawsey0220/anderson/scripts/samtools_array.sbatch .
# modify the array script to match how many samples to map (69, so 0-68)
sbatch samtools_array.sbatch -f ../map_samples.txt
```

Concatenate the summary files together for input to a spreadsheet
```bash
index=1
for sample in $(cat ../map_samples.txt); do
if [ $index == 1 ]; then
head -n 1 "$sample"/"$sample"_samtools_summary.txt > samtools_overall.txt
fi
tail -n 1 "$sample"/"$sample"_samtools_summary.txt >> samtools_overall.txt
index=$((index + 1))
done
# NOTE: the positions of hets and Ns will be extensive, so could drop those columns for limiting file size
# cut -f 8,10 --complement samtools_overall.txt > samtools_overall_smaller.txt
```

Split out the fastas from the consensus sequences and move to each sample directory
```bash
# start an interactive job
salloc -p work -n 1 -N 1 -c 4 -A pawsey0220 --time=01:00:00
module load python/3.11.6

# grab the consensus sequences and move them
for sample in $(cat ../map_samples.txt); do
python3 /software/projects/pawsey0220/anderson/scripts/fasta_splitting.py "$sample"/"$sample"_consensus.fasta
for outfasta in "$sample"-*.fasta; do
mv "$outfasta" ../"$sample"/"${outfasta/${sample}-/consens_}"
done
done

# concatenate across samples per locus, removing locus info from headers
cd ..
mkdir aligning
for locus in $(cat loci.txt); do
for sample in $(cat map_samples.txt); do
if [ -f "$sample"/consens_"$locus".fasta ]; then
sed "/^>/s/-$locus//" "$sample"/consens_"$locus".fasta >> aligning/"$locus".fasta
fi
done
done
exit
```

Run alignments of the loci with MAFFT
```bash
cp /software/projects/pawsey0220/anderson/scripts/align_phylo.sbatch .
# modify the script for the cleaning step to remove positions with > 80% gaps and retain default of sample removal when missing > 75%
sbatch align_phylo.sbatch -a "none" -c "y" -f "aligning" -r "y"
```

Generate a distance matrix from the alignments for use with SplitsTree4
```bash
Rscript ~/scripts/align_to_distance.R -p "g" in_align/*.fasta
```
The resulting distance network is more resolved than from the HybPiper dataset, though it still lacks resoultion outside of SER and CUN  
A few samples suggest errors in assembly and/or mapping: FOR3-01,FOR1-07,FOR1-02  
HERB2 is recovered in a hybrid position (shared distances with CUN), but HERB1 is not  

Note: essentially the same structure is recovered by generating a distance matrix from the exon alignments  
This suggests adding read mapping to full loci is not effectively changing relationships from unambiguous assemblies  
In other words, underlying heterozygosity is not effective at clustering the hybrids for these loci with skimming data  

Heterozygosity is nonetheless clearly different (outliers) for the HERB samples:  
HERB samples: 1.27% heterozygosity (1.25–1.29)  
other samples: 0.41% heterozygosity (0.11–0.64)  
