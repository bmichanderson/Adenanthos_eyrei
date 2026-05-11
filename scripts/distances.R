
##########
# Author: B.M. Anderson
# Date: Nov 2021
# Modified: Oct 2023 (removed plotting and simplified);
#	May 2026 (averaging distances when more than one SNP per locus)
# Description: calculate and output distances between samples from a VCF file
##########


# load required libraries
suppressMessages(library(adegenet))
suppressMessages(library(ape))
suppressMessages(library(pofadinr))
suppressMessages(library(vcfR))


# Define a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate distances from an input VCF file, then output the matrix\n")
		cat("Distances will be average per locus if there are >1 SNPs per locus\n")
		cat("Usage: Rscript distances.R -d dist_method -o output -v vcf_file [-s samples]\n")
		cat("Options:\n")
		cat("\t-d\tThe distance method: \"G\"ENPOFAD [default], \"E\"uclidean, \"M\"ATCHSTATES\n")
		cat("\t-o\tThe output file name prefix [default \"output\"]\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t-s\tA file with sampleID and sample name (tab separated, one per line) [optional]\n")
	} else {
		cat(help_message)
	}
}


# Parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	dist_method <- "G"
	output <- "output"
	vcf_present <- FALSE
	samples_present <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-d") {
			dist_method <- args[index + 1]
		} else if (args[index] == "-o") {
			output <- args[index + 1]
		} else if (args[index] == "-v") {
			vcf_present <- TRUE
			vcf_file <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

if (! vcf_present) {
	stop(help("Missing argument for vcf file!\n"), call. = FALSE)
}


# read in the input
vcf <- read.vcfR(vcf_file, verbose = FALSE)
cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
	length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")

if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
}


# check and record samples present
samples <- colnames(vcf@gt)[2: length(colnames(vcf@gt))]
if (samples_present) {
	if (length(samples) > nrow(sample_table)) {
		cat("Not enough samples in the samples file, so some will not be renamed!\n")
	}
}


# break up into vcfR objects per locus
vcfr_list <- vector("list")
index <- 1
loci <- unique(vcf@fix[, 1])
for (locus in loci) {
	this_vcfr <- vcf[vcf@fix[, 1] == locus, ]
	vcfr_list[[index]] <- this_vcfr
	index <- index + 1
}


# calculate distances between samples per locus
# add the matrix to a list, to be averaged afterwards
dist_list <- vector("list")
index <- 1
for (locus_vcfr in vcfr_list) {
	# calculate the distances
	if (dist_method == "G" || dist_method == "M") {		# GENPOFAD or MATCHSTATES
		dnabin <- vcfR2DNAbin(locus_vcfr, consensus = TRUE, extract.haps = FALSE)
		# change the default behaviour of converting missing into "n"
		temp <- as.character(dnabin)
		temp[temp == "n"] <- "?"
		dnabin <- as.DNAbin(temp)
		if (dist_method == "G") {
			distance <- dist.snp(dnabin, model = "GENPOFAD")
		} else {
			distance <- dist.snp(dnabin, model = "MATCHSTATES")
		}
	} else if (dist_method == "E") {		## Euclidean
		genl <- vcfR2genlight(locus_vcfr)
		distance <- dist(as.matrix(genl))
	} else {
		stop(help("Distance method specified incorrectly!\n"), call. = FALSE)
	}

	# adjust the matrix to ensure all samples are present
	# if there are missing samples, create columns and rows of NAs for them
	dist_mat <- as.matrix(distance)
	for (sample in samples) {
		if (! sample %in% rownames(dist_mat)) {
			curr_names <- rownames(dist_mat)
			dist_mat <- rbind(cbind(dist_mat, rep(NA, nrow(dist_mat))), rep(NA, ncol(dist_mat) + 1))
			rownames(dist_mat) <- c(curr_names, sample)
			colnames(dist_mat) <- c(curr_names, sample)
		}
	}

	# reorder to the order of samples (ensuring same order across all matrices)
	dist_mat <- dist_mat[samples, samples]

	# add to the list
	dist_list[[index]] <- dist_mat

	# increment
	index <- index + 1
}


# calculate the average distance matrix
# since the matrices have the same dimensions and order, calculate means across matrices
# see: https://stackoverflow.com/a/19220503
cat(paste0("\nTransforming into an array..."))
mat_array <- array(unlist(dist_list), c(length(samples), length(samples), length(dist_list)))
cat(paste0("\nCalculating the average distance matrix for ", length(samples), " samples...\n"))
distances <- as.matrix(rowMeans(mat_array, dims = 2, na.rm = TRUE))
rownames(distances) <- samples


# replace sample names if provided
if (samples_present) {
	for (index in seq_len(length(samples))) {
		id <- samples[index]
		if (id %in% sample_table$V1) {
			new_name <- sample_table$V2[match(id, sample_table$V1)]
			rownames(distances)[index] <- new_name
			samples[index] <- new_name
		} else {
			cat(paste0("Sample ", id, " not found in sample table\n"))
		}
	}
}


# if the rownames have spaces, need to quote for Nexus output
space_present <- FALSE
for (rowname in rownames(distances)) {
	if (length(strsplit(rowname, " ")[[1]]) > 1) {
		space_present <- TRUE
		break
	}
}
if (space_present) {
	rownames(distances) <- paste0("'", rownames(distances), "'")
	samples <- paste0("'", samples, "'")
}


# Output the distance matrix
if (dist_method == "G") {
	suffix <- "_GENPOFAD.nex"
} else if (dist_method == "M") {
	suffix <- "_MATCHSTATES.nex"
} else {
	suffix <- "_Euclidean.nex"
}

taxa_block <- paste0("BEGIN TAXA;\n\tDIMENSIONS NTAX=", length(samples), ";\n\t",
	"TAXLABELS ", paste(samples, collapse = " "), ";\nEND;\n")
dist_block <- paste0("BEGIN DISTANCES;\n\tFORMAT\n\t\tTRIANGLE=BOTH\n\t\tDIAGONAL\n\t\t",
	"LABELS=LEFT\n\t;\n\tMATRIX\n")
outfile <- file(paste0(output, suffix), open = "w")
writeLines("#NEXUS", con = outfile)
writeLines(taxa_block, con = outfile)
writeLines(dist_block, con = outfile)
write.table(distances, file = outfile, col.names = FALSE,
	append = TRUE, quote = FALSE)
writeLines("\t;\nEND;\n", con = outfile)
close(outfile)
cat("\n\nDone!\n")
