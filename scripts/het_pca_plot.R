
##########
# Author: B.M. Anderson
# Date: May 2026
# Description: determine heterozygosity and run a PCA for samples from a VCF file
#	Samples can be selected using an optional input samples file (default run on all samples)
#	The samples file can have a second column for point colour (defaul black)
##########


# load required libraries
suppressMessages(library(adegenet))
suppressMessages(library(vcfR))


# Define a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to determine heterozygosity and run a PCA for samples from an input VCF file\n")
		cat("Usage: Rscript het_pca_plot.R [-s samples_file] [-l] vcf_file\n")
		cat("Options:\n")
		cat("\t-l\tFlag to plot labels above the points [default do not]\n")
		cat("\t-o\tThe output file name prefix [default \"output\"]\n")
		cat("\t-s\tA file with sampleID and [optionally] point colour (tab separated, one per line) [optional]\n")
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
	extra <- 1
	catch <- TRUE
	output <- "output"
	samples_present <- FALSE
	plot_labels <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-l") {
			plot_labels <- TRUE
		} else if (args[index] == "-o") {
			output <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
			catch <- FALSE
		} else {
			if (catch) {
				catch_args[extra] <- args[index]
				extra <- extra + 1
			} else {
				catch <- TRUE
			}
		}
	}
}

if (length(catch_args) == 0) {
	stop(help("Missing VCF file!\n"), call. = FALSE)
}


# read in the input (first extra argument, assumed to be the VCF file)
vcf_file <- catch_args[[1]]
vcfr <- read.vcfR(vcf_file, verbose = FALSE)
cat("Read in a VCF with", ncol(vcfr@gt) - 1, "samples,",
	length(unique(vcfr@fix[, 1])), "loci and", nrow(vcfr@fix), "SNPs\n")

samples <- colnames(vcfr@gt)[2: length(colnames(vcfr@gt))]

if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
	if (ncol(sample_table) == 2) {
		colours <- sample_table[, 2]
	} else {
		colours <- rep("black", ncol(sample_table))
	}
} else {
	colours <- rep("black", length(samples))
}


# filter samples if a samples file was provided
if (samples_present) {
	keep_samples <- sample_table[, 1][sample_table[, 1] %in% samples]
	keep_colours <- colours[sample_table[, 1] %in% samples]
	filt_vcfr <- vcfr[sample = keep_samples]
	samples <- keep_samples
	colours <- keep_colours
	vcfr <- filt_vcfr
}


# convert the vcfr to a genlight
genl <- suppressWarnings(vcfR2genlight(vcfr))
if (samples_present) {
	cat("Filtered and ")
}
cat(paste0("Converted the VCF to a genlight with ", nInd(genl), " samples, ",
	length(unique(genl@chromosome)), " loci, and ", nLoc(genl), " SNPs\n"))


# measure heterozygosity for each sample
comp_mat <- as.matrix(genl)
het_df <- as.data.frame(100 * rowSums(comp_mat == 1, na.rm = TRUE) / (ncol(comp_mat) - rowSums(is.na(comp_mat))))
het_df <- as.data.frame(cbind(rownames(het_df), het_df[, 1]))
het_df[, 2] <- as.numeric(het_df[, 2])


# run a PCA
pca <- glPca(genl, nf = 5)
pca_var <- round(100 * pca$eig / sum(pca$eig), digits = 1)


# plot
pdf(paste0(output, "_het_pca.pdf"), width = 7, height = 7)
plot(pca$scores[, 1], het_df[, 2],
	ylim = (c(0, 1.2 * max(het_df[, 2]))),
	xlab = paste0("PC1 (", pca_var[1], "% of variation)"),
	ylab = "SNP % heterozygosity",
	col = colours,
	pch = 19,
	cex = 2
	)

# if desired, plot point labels
if (plot_labels) {
	text(pca$scores[, 1], het_df[, 2], het_df[, 1], pos = 3)
}

invisible(dev.off())
